#!/usr/bin/env bash
# Global Ancestry Inference pipeline using the gnomAD HGDP+1kGP subset
# (genomes only) as the reference panel.
#
# Reference panel: ~4,094 individuals from HGDP + 1000 Genomes Project,
# re-sequenced to high coverage and re-called through the gnomAD pipeline.
# ~80 sub-populations with ethnographic labels (YRI, CEU, Sardinian, Karitiana, ...).
#
# WHY v3.1.2 paths even though we want "v4.x semantics":
#   gnomAD v4.0 (Nov 2023) and v4.1 (April 2024) added a new exome callset and
#   refreshed site-level allele-frequency annotations, but they did NOT
#   republish the HGDP+1kGP per-sample dense VCFs — those individual-level
#   genotypes remain at the v3.1.2 release path. Verified empirically:
#     gsutil ls gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/ | grep hgdp
#       (returns nothing)
#     gsutil ls gs://gcp-public-data--gnomad/release/4.0/vcf/genomes/ | grep hgdp
#       (returns nothing)
#     gsutil ls gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/ | grep hgdp
#       (returns the 22 autosomes + X/Y + sample-meta TSV)
#   PCA / ADMIXTURE require per-sample genotypes, so we use the v3.1.2 path.
#   This is the same genotype data that v4.x conceptually inherits — not a
#   downgrade.
#
# Public location:
#   gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/
#     gnomad.genomes.v3.1.2.hgdp_tgp.chr{1..22}.vcf.bgz
#     gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz
#   (HTTPS mirror via storage.googleapis.com)
#
# IMPORTANT:
#   - Per-chromosome dense VCFs are large (multi-GB compressed). We stream with
#     `bcftools view -R` so only the sites you need are transferred.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="${PROJECT_DIR}"
PIPELINE_DIR="${SCRIPT_DIR}"
SCRIPTS_DIR="${PIPELINE_DIR}/scripts"
WORKING="${PIPELINE_DIR}/working"
RESULTS="${PIPELINE_DIR}/results"
REFDIR="${PIPELINE_DIR}/reference"
TOOLS_DIR="${PIPELINE_DIR}/tools"
CONDA_ENV="ancestry_pipeline"
THREADS=$(sysctl -n hw.ncpu || nproc || echo 4)

# --- QC thresholds (same as 1KG pipeline) ---
GENO=0.05
MIND=0.1
MAF=0.01
HWE=1e-4

# --- LD pruning ---
LD_WINDOW=50
LD_STEP=5
LD_R2=0.2

# --- ADMIXTURE K range ---
# HGDP+1kGP resolves more structure, so we extend the upper K.
# Typical "interesting" Ks for this panel: 5–10.
K_MIN=3
K_MAX=10
MIN_GS=0.15

# ============================================================================
# gnomAD HGDP+1kGP subset paths
# ============================================================================
# v3.1.2 is the only release that has per-sample dense VCFs for HGDP+1kGP.
# v4.0 and v4.1 do not republish them (verified — see top-of-file note).
# Verify with:
#   gsutil ls gs://gcp-public-data--gnomad/release/${GNOMAD_VERSION}/vcf/genomes/ | grep hgdp
GNOMAD_VERSION="3.1.2"

# Mirror selection: AWS S3 mirror is HTTP/1.1, GCS mirror is HTTP/2.
# The HTTP/2 framing layer (libcurl error 16) corrupts byte-range reads at
# multi-GB offsets and bcftools silently truncates output — verified with a
# 30-region pull where GCS returned 15 KB vs AWS 603 KB for the same input.
# Both mirrors host the identical files (same ETag suffix, same release).
GNOMAD_BASE="https://gnomad-public-us-east-1.s3.amazonaws.com/release/${GNOMAD_VERSION}/vcf/genomes"
GNOMAD_VCF_TPL="${GNOMAD_BASE}/gnomad.genomes.v${GNOMAD_VERSION}.hgdp_tgp.chr{CHR}.vcf.bgz"

# Sample metadata (population labels). Co-located with the per-chr VCFs in the
# bucket. File is a bgzipped TSV — we decompress it after download.
# Verified for v3.1.2 path:
#   gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz
GNOMAD_META_URL="${GNOMAD_BASE}/gnomad.genomes.v${GNOMAD_VERSION}.hgdp_1kg_subset_sample_meta.tsv.bgz"

log()  { echo "[$(date '+%H:%M:%S')] $*"; }
warn() { echo "[$(date '+%H:%M:%S')] WARNING: $*" >&2; }
die()  { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }
step_done() { [ -f "${WORKING}/$1.done" ]; }
mark_done() { touch "${WORKING}/$1.done"; }

# ============================================================================
# Step 0: setup
# ============================================================================
step0_setup() {
    if step_done "step0"; then
        log "Step 0 already done, skipping setup"
        return
    fi

    log "STEP 0: Setup"
    mkdir -p "${WORKING}" "${RESULTS}/pca" "${RESULTS}/admixture" "${RESULTS}/plots" \
             "${REFDIR}" "${TOOLS_DIR}"

    CONDA_EXE="$(which conda || echo "")"
    if [ -z "$CONDA_EXE" ]; then
        die "conda not found. Please install miniconda/miniforge first."
    fi
    CONDA_BASE="$(dirname "$(dirname "$CONDA_EXE")")"
    source "${CONDA_BASE}/etc/profile.d/conda.sh" || true

    if conda info --envs | grep -q "${CONDA_ENV}"; then
        log "Conda env exists: ${CONDA_ENV}"
    else
        log "Creating conda env: ${CONDA_ENV}"
        conda create -n "${CONDA_ENV}" -c bioconda -c conda-forge \
            python=3.10 plink2 bcftools samtools htslib \
            matplotlib pandas numpy scipy admixture \
            -y
    fi
    conda activate "${CONDA_ENV}"

    for tool in plink2 bcftools samtools python3 curl; do
        command -v "$tool" >/dev/null 2>&1 || die "$tool not found"
    done

    if ! command -v admixture >/dev/null 2>&1; then
        warn "ADMIXTURE not found from conda. PCA will still run, but ADMIXTURE may be skipped."
    fi

    mark_done "step0"
}

# ============================================================================
# Step 1: convert DDNA -> PLINK
# ============================================================================
step1_convert() {
    if step_done "step1"; then
        log "Step 1 already done, skipping conversion"
        return
    fi

    log "STEP 1: Convert DDNA files to PLINK"

    python3 "${SCRIPTS_DIR}/convert_ddna_to_plink.py" \
        "${DATA_DIR}" "${WORKING}/study_raw" "${MIN_GS}"

    # Drop multiallelic SNPs from the DDNA-derived TPED before plink2 conversion.
    python3 - <<FILTER_TPED
from pathlib import Path

working = Path("${WORKING}")
inp = working / "study_raw.tped"
out = working / "study_raw.biallelic.tped"
bad = working / "multiallelic_snps_removed.txt"

n_total = n_keep = n_bad = 0
with inp.open() as f, out.open("w") as fo, bad.open("w") as fb:
    for line in f:
        n_total += 1
        parts = line.rstrip("\n").split()
        snp = parts[1]
        alleles = parts[4:]
        observed = {a for a in alleles if a not in {"0", "N", "-", "."}}
        if len(observed) <= 2:
            fo.write(line)
            n_keep += 1
        else:
            fb.write(f"{snp}\t{','.join(sorted(observed))}\n")
            n_bad += 1

print(f"Filtered TPED: kept {n_keep}/{n_total} biallelic SNPs; removed {n_bad} multiallelic SNPs")
FILTER_TPED

    plink2 --tped "${WORKING}/study_raw.biallelic.tped" \
           --tfam "${WORKING}/study_raw.tfam" \
           --make-bed \
           --out "${WORKING}/study_raw" \
           --threads "${THREADS}" \
           --allow-extra-chr

    rm -f "${WORKING}/study_raw.tped"

    log "Converted: $(wc -l < "${WORKING}/study_raw.bim") SNPs, $(wc -l < "${WORKING}/study_raw.fam") samples"
    mark_done "step1"
}

# ============================================================================
# Step 2: QC
# ============================================================================
step2_qc() {
    if step_done "step2"; then
        log "Step 2 already done, skipping QC"
        return
    fi

    log "STEP 2: QC"

    plink2 --bfile "${WORKING}/study_raw" \
           --rm-dup exclude-all \
           --make-bed \
           --out "${WORKING}/study_nodup" \
           --threads "${THREADS}"

    plink2 --bfile "${WORKING}/study_nodup" \
           --geno "${GENO}" \
           --mind "${MIND}" \
           --maf "${MAF}" \
           --hwe "${HWE}" \
           --make-bed \
           --out "${WORKING}/study_qc" \
           --threads "${THREADS}"

    {
        echo "=== QC Report ==="
        echo "Input SNPs: $(wc -l < "${WORKING}/study_raw.bim")"
        echo "Input samples: $(wc -l < "${WORKING}/study_raw.fam")"
        echo "Final SNPs: $(wc -l < "${WORKING}/study_qc.bim")"
        echo "Final samples: $(wc -l < "${WORKING}/study_qc.fam")"
        echo "geno=${GENO}, mind=${MIND}, maf=${MAF}, hwe=${HWE}"
    } > "${RESULTS}/qc_report.txt"

    rm -f "${WORKING}"/study_nodup.*
    mark_done "step2"
}

# ============================================================================
# Step 2b: shrink the study BIM to gnomAD's pre-pruned ancestry-informative SNPs
# ============================================================================
# Why this matters: Step 3 sends one BED region per study SNP to the remote
# gnomAD VCFs. The post-QC array has ~400-700k SNPs; intersecting with
# gnomAD's PCA loadings shrinks that to ~30-80k variants — ~5-15x less data
# to stream and ~5-15x less wall time.
#
# Why not local LD pruning (--indep-pairwise)? PLINK refuses with N<50, and
# rightfully so: r² estimated from 10 samples is noise, so the "pruned" set
# would be a random subset. The loadings HT IS the LD-pruned set computed
# by gnomAD on all 4094 reference samples, AND it's been filtered for
# ancestry-informativeness across global populations. Strictly better than
# anything we could compute locally.
#
# Step 5 still re-prunes on the joint (study∩ref) merge as a safety net.
LOADINGS_HT="${LOADINGS_HT:-${REFDIR}/pca_loadings/gnomad.v3.1.pca_loadings.ht}"

step2b_ld_prune_study() {
    if step_done "step2b"; then
        log "Step 2b already done, skipping loadings-intersect"
        return
    fi

    log "STEP 2b: Intersect study BIM with gnomAD PCA loadings variants"

    LOADINGS_TSV="${WORKING}/loadings_variants.tsv"
    if [ ! -s "${LOADINGS_TSV}" ]; then
        if [ ! -d "${LOADINGS_HT}" ]; then
            warn "Loadings HT not found: ${LOADINGS_HT}"
            warn "Falling back to un-pruned study_qc.bim for download set."
            warn "(This will be slow. Run Path A first or set LOADINGS_HT.)"
            for ext in bed bim fam; do
                cp "${WORKING}/study_qc.${ext}" "${WORKING}/study_qc_pruned.${ext}"
            done
            mark_done "step2b"
            return
        fi
        log "Extracting variant list from ${LOADINGS_HT}"
        python "${SCRIPTS_DIR}/extract_loadings_variants.py" \
               "${LOADINGS_HT}" "${LOADINGS_TSV}"
    fi
    log "Loadings variants available: $(wc -l < "${LOADINGS_TSV}")"

    # Rename study variants to chrom:pos (matches col 5 of loadings TSV).
    plink2 --bfile "${WORKING}/study_qc" \
           --set-all-var-ids '@:#' \
           --new-id-max-allele-len 10 truncate \
           --make-bed \
           --out "${WORKING}/study_qc_chrpos_for_prune" \
           --threads "${THREADS}"

    # Build the extract list: variant IDs that appear in BOTH the study BIM
    # and the loadings variant set.
    awk '{print $5}' "${LOADINGS_TSV}" | sort -u > "${WORKING}/loadings_ids.txt"
    awk '{print $2}' "${WORKING}/study_qc_chrpos_for_prune.bim" | sort -u \
        > "${WORKING}/study_chrpos_ids.txt"
    comm -12 "${WORKING}/loadings_ids.txt" "${WORKING}/study_chrpos_ids.txt" \
        > "${WORKING}/study_loadings_intersect_ids.txt"

    N_INT=$(wc -l < "${WORKING}/study_loadings_intersect_ids.txt")
    [ "${N_INT}" -ge 1000 ] || die "Only ${N_INT} variants intersect with loadings — check builds (study should be GRCh38)."

    plink2 --bfile "${WORKING}/study_qc_chrpos_for_prune" \
           --extract "${WORKING}/study_loadings_intersect_ids.txt" \
           --make-bed \
           --out "${WORKING}/study_qc_pruned" \
           --threads "${THREADS}"

    N_BEFORE=$(wc -l < "${WORKING}/study_qc.bim")
    N_AFTER=$(wc -l < "${WORKING}/study_qc_pruned.bim")
    log "Study SNPs: ${N_BEFORE} -> ${N_AFTER} after loadings-intersect (~$(( 100 * N_AFTER / N_BEFORE ))%)"

    rm -f "${WORKING}"/study_qc_chrpos_for_prune.* \
          "${WORKING}/loadings_ids.txt" "${WORKING}/study_chrpos_ids.txt"
    mark_done "step2b"
}

# ============================================================================
# Step 3: download gnomAD HGDP+1kGP reference (sample meta + per-chr VCF slices)
# ============================================================================
step3_download_reference() {
    if step_done "step3"; then
        log "Step 3 already done, skipping reference download"
        return
    fi

    log "STEP 3: Download gnomAD v${GNOMAD_VERSION} HGDP+1kGP reference"

    # --- Sample metadata (bgzipped TSV co-located with the per-chr VCFs) ---
    META_TSV="${REFDIR}/hgdp_tgp_sample_meta.tsv"
    META_BGZ="${META_TSV}.bgz"
    if [ ! -s "${META_TSV}" ]; then
        log "Downloading sample metadata"
        if ! curl -fL "${GNOMAD_META_URL}" -o "${META_BGZ}.tmp"; then
            warn "Could not download ${GNOMAD_META_URL}"
            warn "List the bucket and update GNOMAD_META_URL at the top of this script:"
            warn "  gsutil ls gs://gcp-public-data--gnomad/release/${GNOMAD_VERSION}/vcf/genomes/ | grep meta"
            warn "Or fall back to v3.1.2 by setting GNOMAD_VERSION=3.1.2"
            die "Sample metadata download failed"
        fi
        mv "${META_BGZ}.tmp" "${META_BGZ}"

        # Decompress bgzipped TSV. bgzip ships with htslib (already in conda env).
        if command -v bgzip >/dev/null 2>&1; then
            bgzip -dc "${META_BGZ}" > "${META_TSV}"
        else
            # Fallback: gzip works too (bgzip is gzip-compatible at the block level)
            gzip -dc "${META_BGZ}" > "${META_TSV}"
        fi
        rm -f "${META_BGZ}"
    fi

    # Build the panel file used by visualization (IID -> pop, super_pop).
    # gnomAD v3.1.2 metadata stores HGDP/1kGP labels inside JSON column: hgdp_tgp_meta
    python3 - <<META_PARSE
import csv, sys, json
from pathlib import Path

inp = Path("${META_TSV}")
out = Path("${REFDIR}/panel_hgdp_tgp.tsv")

n = 0
n_skip = 0

with inp.open() as f, out.open("w") as fo:
    reader = csv.DictReader(f, delimiter="\t")
    fields = reader.fieldnames or []

    if "s" not in fields:
        sys.stderr.write(f"Could not find sample ID column 's'. Saw: {fields[:20]}\n")
        sys.exit(1)

    if "hgdp_tgp_meta" not in fields:
        sys.stderr.write(f"Could not find JSON column 'hgdp_tgp_meta'. Saw: {fields[:20]}\n")
        sys.exit(1)

    fo.write("sample\tpop\tsuper_pop\tproject\n")

    for row in reader:
        sid = (row.get("s") or "").strip()
        meta_raw = (row.get("hgdp_tgp_meta") or "").strip()

        if not sid or not meta_raw:
            n_skip += 1
            continue

        try:
            meta = json.loads(meta_raw)
        except Exception:
            n_skip += 1
            continue

        pop = meta.get("population")
        sup = meta.get("genetic_region")
        proj = meta.get("project")

        pop = "" if pop is None else str(pop).strip()
        sup = "" if sup is None else str(sup).strip()
        proj = "" if proj is None else str(proj).strip()

        # Skip synthetic / unlabeled samples
        if not pop or pop.lower() in {"none", "null", "na"}:
            n_skip += 1
            continue

        # If genetic_region is missing, use project as fallback group label
        if not sup or sup.lower() in {"none", "null", "na"}:
            sup = proj if proj else "unknown"

        fo.write(f"{sid}\t{pop}\t{sup}\t{proj}\n")
        n += 1

print(f"Wrote {n} labeled samples to {out}; skipped {n_skip} unlabeled/malformed samples")
META_PARSE

    # Sample list for bcftools -S (just the IDs)
    awk 'NR>1 {print $1}' "${REFDIR}/panel_hgdp_tgp.tsv" > "${REFDIR}/hgdp_tgp_samples.txt"
    log "Reference samples: $(wc -l < "${REFDIR}/hgdp_tgp_samples.txt")"

    # --- Build a BED of study SNP positions to slice gnomAD by ---
    # Use the LD-pruned BIM if present (~3-6x fewer rows than study_qc.bim).
    POSITIONS_SRC="${WORKING}/study_qc_pruned.bim"
    if [ ! -s "${POSITIONS_SRC}" ]; then
        warn "study_qc_pruned.bim not found; falling back to un-pruned study_qc.bim"
        POSITIONS_SRC="${WORKING}/study_qc.bim"
    fi
    log "BED-positions source: ${POSITIONS_SRC} ($(wc -l < "${POSITIONS_SRC}") SNPs)"

    awk '{print "chr"$1"\t"$4-1"\t"$4"\t"$2}' "${POSITIONS_SRC}" \
        | sort -k1,1 -k2,2n > "${WORKING}/study_positions_chr.bed"

    # Cleanup .tmp files from interrupted prior runs.
    find "${REFDIR}" -name 'hgdp_tgp_study_snps_chr*.vcf.gz.tmp*' -delete 2>/dev/null || true

    # --- Pre-download tabix index files (.tbi) for all 22 chromosomes -------
    # Without local indexes, bcftools falls back to LINEAR streaming of the
    # remote VCF (gigabytes per chr) which inevitably blows the HTTP/2
    # connection mid-stream (libcurl error 16 / framing layer). With a local
    # .tbi, bcftools does precise byte-range fetches on the remote VCF — each
    # request is tiny and short-lived, completely avoiding the long-stream
    # issue. Each .tbi is ~30-100 KB; total ~1 MB.
    TBI_DIR="${REFDIR}/tbi"
    mkdir -p "${TBI_DIR}"
    log "Pre-downloading .tbi indexes for all 22 chrs (one-time, <1 MB total)"
    for CHR in $(seq 1 22); do
        local TBI_LOCAL="${TBI_DIR}/gnomad.genomes.v${GNOMAD_VERSION}.hgdp_tgp.chr${CHR}.vcf.bgz.tbi"
        local TBI_URL="${GNOMAD_VCF_TPL/\{CHR\}/${CHR}}.tbi"
        if [ ! -s "${TBI_LOCAL}" ]; then
            for tbi_try in 1 2 3; do
                if curl -fsSL --http1.1 -o "${TBI_LOCAL}.tmp" "${TBI_URL}"; then
                    mv "${TBI_LOCAL}.tmp" "${TBI_LOCAL}"
                    break
                fi
                rm -f "${TBI_LOCAL}.tmp"
                sleep 5
            done
            [ -s "${TBI_LOCAL}" ] || die "Failed to download .tbi for chr${CHR}"
        fi
    done
    log "  .tbi indexes ready: $(ls "${TBI_DIR}" | wc -l | tr -d ' ') files"

    # Number of chromosomes to download in parallel. With unstable HTTPS,
    # 1-2 is the sweet spot — gnomAD/GCS throttles per-client concurrent range
    # requests, and a single failed connection takes out the whole job.
    PARALLEL_CHRS="${PARALLEL_CHRS:-2}"

    # Coalesce 1bp BED regions within MERGE_DIST bp into single regions.
    # 30k 1bp regions -> ~5k merged regions, ~6-10x fewer HTTP requests.
    MERGE_DIST="${MERGE_DIST:-10000}"

    # Retry the remote bcftools fetch this many times before giving up
    # (with linear backoff: 30s, 60s, 90s, ...).
    MAX_TRIES="${MAX_TRIES:-5}"

    # gnomAD's HTTPS endpoint frequently drops HTTP/2 connections mid-stream
    # (libcurl error 16 "HTTP2 framing layer"). Splitting each chr into many
    # small bcftools requests means a single connection drop only loses a
    # chunk, not the whole chr. 30-50 regions per chunk is a good balance:
    # small enough to almost always complete, large enough that TLS handshake
    # overhead doesn't dominate.
    CHUNK_REGIONS="${CHUNK_REGIONS:-30}"

    download_one_chr() {
        local CHR=$1
        local CHR_VCF="${REFDIR}/hgdp_tgp_study_snps_chr${CHR}.vcf.gz"
        local POS_BED="${WORKING}/positions_chr${CHR}.bed"
        local POS_MERGED="${WORKING}/positions_chr${CHR}.merged.bed"

        if [ -s "${CHR_VCF}" ]; then
            log "chr${CHR}: already exists"
            return 0
        fi

        awk -v chr="chr${CHR}" '$1 == chr' "${WORKING}/study_positions_chr.bed" > "${POS_BED}"
        local N_POS
        N_POS=$(wc -l < "${POS_BED}")
        if [ "${N_POS}" -eq 0 ]; then
            log "chr${CHR}: no study SNPs"
            return 0
        fi

        # Coalesce nearby SNPs to reduce HTTP round-trips. The trade-off:
        # each request fetches a slightly larger contiguous region, but we
        # pay HTTP+TLS overhead far fewer times.
        sort -k1,1 -k2,2n "${POS_BED}" | awk -v dist="${MERGE_DIST}" -v OFS='\t' '
            NR == 1 { chr=$1; s=$2; e=$3; next }
            $1 == chr && $2 - e <= dist { if ($3 > e) e = $3; next }
            { print chr, s, e; chr=$1; s=$2; e=$3 }
            END { if (chr != "") print chr, s, e }
        ' > "${POS_MERGED}"

        local N_MERGED
        N_MERGED=$(wc -l < "${POS_MERGED}")

        # Use ##idx## to point bcftools at the LOCAL .tbi for this chr.
        # Without this, bcftools fails to load remote .tbi via HTTP/2 and
        # silently falls back to linear streaming the entire 30+ GB VCF —
        # which is what was causing the "framing layer" errors at huge offsets.
        local TBI_LOCAL="${TBI_DIR}/gnomad.genomes.v${GNOMAD_VERSION}.hgdp_tgp.chr${CHR}.vcf.bgz.tbi"
        local SRC_VCF="${GNOMAD_VCF_TPL/\{CHR\}/${CHR}}##idx##${TBI_LOCAL}"

        # Split POS_MERGED into chunks of CHUNK_REGIONS regions each. Each
        # chunk becomes its own bcftools fetch — short-lived requests survive
        # gnomAD's flaky HTTP/2.
        local CHUNK_DIR="${WORKING}/chunks_chr${CHR}"
        mkdir -p "${CHUNK_DIR}"
        # Don't wipe — preserve any chunks downloaded by a prior interrupted run.
        # split overwrites the region-list files (regions_aaaa, etc.) deterministically
        # with the same content, so old + new region lists are identical.
        split -a 4 -l "${CHUNK_REGIONS}" "${POS_MERGED}" "${CHUNK_DIR}/regions_"
        local N_CHUNKS
        # Glob `regions_????` matches the 4-char split suffix only, excluding
        # any `regions_aaaa.vcf.gz` produced by prior runs.
        N_CHUNKS=$(ls "${CHUNK_DIR}"/regions_???? 2>/dev/null | wc -l | tr -d ' ')
        log "chr${CHR}: ${N_POS} SNPs -> ${N_MERGED} regions -> ${N_CHUNKS} chunks (${CHUNK_REGIONS}/chunk)"

        local CONCAT_LIST="${CHUNK_DIR}/concat.list"
        > "${CONCAT_LIST}"
        local CHUNK_OK=0
        local CHUNK_FAIL=0

        for CHUNK in "${CHUNK_DIR}"/regions_????; do
            local CHUNK_VCF="${CHUNK}.vcf.gz"
            # Resume: a non-empty .vcf.gz from a prior run is already valid output.
            if [ -s "${CHUNK_VCF}" ]; then
                echo "${CHUNK_VCF}" >> "${CONCAT_LIST}"
                CHUNK_OK=$((CHUNK_OK + 1))
                continue
            fi
            local TRY=1
            local CHUNK_DONE=0
            while [ "${TRY}" -le "${MAX_TRIES}" ]; do
                if bcftools view \
                        -R "${CHUNK}" \
                        -f PASS \
                        -v snps -m 2 -M 2 \
                        "${SRC_VCF}" 2> "${CHUNK_VCF}.err" \
                    | bcftools annotate \
                        -x 'INFO,^FORMAT/GT' \
                        -Oz -o "${CHUNK_VCF}" 2>> "${CHUNK_VCF}.err"; then
                    CHUNK_DONE=1
                    rm -f "${CHUNK_VCF}.err"
                    break
                fi
                local DELAY=$((TRY * 15))
                warn "chr${CHR} $(basename "${CHUNK}"): attempt ${TRY}/${MAX_TRIES} failed; retry in ${DELAY}s"
                tail -n 2 "${CHUNK_VCF}.err" >&2 2>/dev/null || true
                rm -f "${CHUNK_VCF}"
                sleep "${DELAY}"
                TRY=$((TRY + 1))
            done

            if [ "${CHUNK_DONE}" -eq 1 ]; then
                echo "${CHUNK_VCF}" >> "${CONCAT_LIST}"
                CHUNK_OK=$((CHUNK_OK + 1))
            else
                warn "chr${CHR} $(basename "${CHUNK}"): all attempts failed; skipping"
                CHUNK_FAIL=$((CHUNK_FAIL + 1))
            fi
        done

        if [ "${CHUNK_OK}" -eq 0 ]; then
            warn "chr${CHR}: no chunks succeeded; chr is empty"
            rm -rf "${CHUNK_DIR}" "${POS_BED}" "${POS_MERGED}"
            return 1
        fi

        log "chr${CHR}: ${CHUNK_OK}/${N_CHUNKS} chunks ok (${CHUNK_FAIL} failed); concatenating"

        # Concat all chunk VCFs into one chr-level VCF.
        # No `-a` because chunks come from a sorted, non-overlapping BED split
        # so no overlap handling is needed — and `-a` would require each chunk
        # be tabix-indexed (we'd have to index 100+ chunks for nothing).
        bcftools concat -f "${CONCAT_LIST}" -Oz -o "${CHR_VCF}.tmp"

        # Filter back down to the exact study positions (in case merge windows
        # pulled in off-target variants). Use `-T` (targets) not `-R` (regions)
        # because the .tmp file is unindexed.
        bcftools view -T "${POS_BED}" "${CHR_VCF}.tmp" -Oz -o "${CHR_VCF}.tmp2"
        mv "${CHR_VCF}.tmp2" "${CHR_VCF}"
        rm -f "${CHR_VCF}.tmp"
        bcftools index -t "${CHR_VCF}" || true
        log "chr${CHR}: done ($(bcftools view -H "${CHR_VCF}" | wc -l | tr -d ' ') sites)"

        rm -rf "${CHUNK_DIR}" "${POS_BED}" "${POS_MERGED}"
    }

    # Spawn downloads with a simple jobs-based throttle (bash 3.2 compatible).
    log "Downloading chrs with PARALLEL_CHRS=${PARALLEL_CHRS}"
    for CHR in $(seq 1 22); do
        while [ "$(jobs -rp | wc -l | tr -d ' ')" -ge "${PARALLEL_CHRS}" ]; do
            sleep 5
        done
        download_one_chr "${CHR}" &
    done
    wait

    mark_done "step3"
}

# ============================================================================
# Step 4: convert per-chr VCFs to a single PLINK fileset
# ============================================================================
step4_prepare_reference() {
    if step_done "step4"; then
        log "Step 4 already done, skipping reference preparation"
        return
    fi

    log "STEP 4: Prepare reference PLINK files"

    MERGE_LIST="${WORKING}/hgdp_merge_list.txt"
    > "${MERGE_LIST}"

    for CHR in $(seq 1 22); do
        CHR_VCF="${REFDIR}/hgdp_tgp_study_snps_chr${CHR}.vcf.gz"
        CHR_BED="${WORKING}/hgdp_chr${CHR}"
        [ -s "${CHR_VCF}" ] || continue

        plink2 --vcf "${CHR_VCF}" \
               --output-chr 26 \
               --set-all-var-ids '@:#' \
               --new-id-max-allele-len 10 truncate \
               --make-bed \
               --out "${CHR_BED}" \
               --threads "${THREADS}"

        echo "${CHR_BED}" >> "${MERGE_LIST}"
    done

    N_CHR=$(wc -l < "${MERGE_LIST}")
    [ "$N_CHR" -gt 0 ] || die "No HGDP+1kGP chromosome files were prepared"

    if [ "$N_CHR" -eq 1 ]; then
        FIRST=$(head -1 "${MERGE_LIST}")
        cp "${FIRST}.bed" "${WORKING}/hgdp_prepared.bed"
        cp "${FIRST}.bim" "${WORKING}/hgdp_prepared.bim"
        cp "${FIRST}.fam" "${WORKING}/hgdp_prepared.fam"
    else
        FIRST=$(head -1 "${MERGE_LIST}")
        REST="${WORKING}/hgdp_merge_rest.txt"
        tail -n +2 "${MERGE_LIST}" > "${REST}"

        plink2 --bfile "${FIRST}" \
               --pmerge-list "${REST}" bfile \
               --make-bed \
               --out "${WORKING}/hgdp_prepared" \
               --threads "${THREADS}"
    fi

    log "Reference: $(wc -l < "${WORKING}/hgdp_prepared.bim") SNPs, $(wc -l < "${WORKING}/hgdp_prepared.fam") samples"

    rm -f "${WORKING}"/hgdp_chr*.bed "${WORKING}"/hgdp_chr*.bim "${WORKING}"/hgdp_chr*.fam
    mark_done "step4"
}

# ============================================================================
# Step 5: merge study + reference, drop ambiguous SNPs, LD-prune
# ============================================================================
step5_merge_prune() {
    if step_done "step5"; then
        log "Step 5 already done, skipping merge/prune"
        return
    fi

    log "STEP 5: Merge study + reference and LD prune"

    plink2 --bfile "${WORKING}/study_qc" \
           --set-all-var-ids '@:#' \
           --make-bed \
           --out "${WORKING}/study_qc_chrpos" \
           --threads "${THREADS}"

    awk '{print $2}' "${WORKING}/study_qc_chrpos.bim" | sort > "${WORKING}/study_snp_ids.txt"
    awk '{print $2}' "${WORKING}/hgdp_prepared.bim"   | sort > "${WORKING}/hgdp_snp_ids.txt"

    comm -12 "${WORKING}/study_snp_ids.txt" "${WORKING}/hgdp_snp_ids.txt" \
        > "${WORKING}/common_snps.txt"

    N_COMMON=$(wc -l < "${WORKING}/common_snps.txt")
    log "Common SNPs (study ∩ HGDP+1kGP): ${N_COMMON}"
    [ "$N_COMMON" -ge 1000 ] || die "Too few common SNPs. Check genome build / chromosome naming."

    plink2 --bfile "${WORKING}/study_qc_chrpos" \
           --extract "${WORKING}/common_snps.txt" \
           --make-bed \
           --out "${WORKING}/study_common" \
           --threads "${THREADS}"

    plink2 --bfile "${WORKING}/hgdp_prepared" \
           --extract "${WORKING}/common_snps.txt" \
           --make-bed \
           --out "${WORKING}/hgdp_common" \
           --threads "${THREADS}"

    # Drop strand-ambiguous (A/T, C/G) SNPs
    awk '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || \
         ($5 == "C" && $6 == "G") || ($5 == "G" && $6 == "C") {print $2}' \
        "${WORKING}/study_common.bim" > "${WORKING}/ambiguous_snps.txt"

    if [ -s "${WORKING}/ambiguous_snps.txt" ]; then
        plink2 --bfile "${WORKING}/study_common" --exclude "${WORKING}/ambiguous_snps.txt" \
               --make-bed --out "${WORKING}/study_noambi" --threads "${THREADS}"
        plink2 --bfile "${WORKING}/hgdp_common"  --exclude "${WORKING}/ambiguous_snps.txt" \
               --make-bed --out "${WORKING}/hgdp_noambi"  --threads "${THREADS}"
    else
        for ext in bed bim fam; do
            cp "${WORKING}/study_common.${ext}" "${WORKING}/study_noambi.${ext}"
            cp "${WORKING}/hgdp_common.${ext}"  "${WORKING}/hgdp_noambi.${ext}"
        done
    fi

    # PLINK 2's --pmerge-list errors with "Non-concatenating --pmerge[-list]
    # is under development" when merging different sample sets (study + ref).
    # That mode of merge is genuinely not implemented in plink2 yet.
    # Use PLINK 1.9's mature --bmerge instead — only needed for this one step.
    # On allele mismatch, plink writes <out>-merge.missnp; we exclude and retry.
    if ! plink --bfile "${WORKING}/study_noambi" \
               --bmerge "${WORKING}/hgdp_noambi" \
               --make-bed \
               --allow-no-sex \
               --out "${WORKING}/merged" \
               --threads "${THREADS}"; then
        warn "Merge failed. Trying to remove mismatched SNPs."
        MISSNP="${WORKING}/merged-merge.missnp"
        if [ -f "$MISSNP" ]; then
            plink --bfile "${WORKING}/study_noambi" --exclude "$MISSNP" \
                  --make-bed --allow-no-sex --out "${WORKING}/study_noambi2" \
                  --threads "${THREADS}"
            plink --bfile "${WORKING}/hgdp_noambi"  --exclude "$MISSNP" \
                  --make-bed --allow-no-sex --out "${WORKING}/hgdp_noambi2" \
                  --threads "${THREADS}"
            plink --bfile "${WORKING}/study_noambi2" \
                  --bmerge "${WORKING}/hgdp_noambi2" \
                  --make-bed --allow-no-sex --out "${WORKING}/merged" \
                  --threads "${THREADS}"
        else
            die "Merge failed and no missnp file was produced."
        fi
    fi

    log "Merged: $(wc -l < "${WORKING}/merged.bim") SNPs, $(wc -l < "${WORKING}/merged.fam") samples"

    plink2 --bfile "${WORKING}/merged" \
           --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
           --out "${WORKING}/merged_ldprune" \
           --threads "${THREADS}"

    plink2 --bfile "${WORKING}/merged" \
           --extract "${WORKING}/merged_ldprune.prune.in" \
           --make-bed \
           --out "${WORKING}/merged_pruned" \
           --threads "${THREADS}"

    log "LD-pruned SNPs: $(wc -l < "${WORKING}/merged_pruned.bim")"
    mark_done "step5"
}

# ============================================================================
# Step 6: PCA
# ============================================================================
step6_pca() {
    if step_done "step6"; then
        log "Step 6 already done, skipping PCA"
        return
    fi

    log "STEP 6: PCA"
    plink2 --bfile "${WORKING}/merged_pruned" \
           --pca 10 \
           --out "${RESULTS}/pca/merged_pca" \
           --threads "${THREADS}"

    mark_done "step6"
}

# ============================================================================
# Step 7: ADMIXTURE
# ============================================================================
step7_admixture() {
    if step_done "step7"; then
        log "Step 7 already done, skipping ADMIXTURE"
        return
    fi

    log "STEP 7: ADMIXTURE"
    if ! command -v admixture >/dev/null 2>&1; then
        warn "ADMIXTURE not available. Skipping."
        mark_done "step7"
        return
    fi

    cp "${WORKING}/merged_pruned.bed" "${RESULTS}/admixture/study_admixture.bed"
    cp "${WORKING}/merged_pruned.bim" "${RESULTS}/admixture/study_admixture.bim"
    cp "${WORKING}/merged_pruned.fam" "${RESULTS}/admixture/study_admixture.fam"

    cd "${RESULTS}/admixture"
    > cv_errors.txt

    for K in $(seq "${K_MIN}" "${K_MAX}"); do
        log "Running ADMIXTURE K=${K}"
        admixture --cv study_admixture.bed "${K}" -j"${THREADS}" 2>&1 | tee "admixture_K${K}.log" || {
            warn "ADMIXTURE K=${K} failed"
            continue
        }
        CV_ERR=$(grep -i "CV error" "admixture_K${K}.log" | awk '{print $NF}' | tail -1)
        [ -n "$CV_ERR" ] && echo "K=${K}: ${CV_ERR}" >> cv_errors.txt
    done

    cd "${PIPELINE_DIR}"
    mark_done "step7"
}

# ============================================================================
# Step 8: visualize
# ============================================================================
step8_visualize() {
    log "STEP 8: Visualization"
    if [ -f "${SCRIPTS_DIR}/visualize_results.py" ]; then
        python3 "${SCRIPTS_DIR}/visualize_results.py" "${RESULTS}" "${REFDIR}" \
            || warn "Visualization failed, but PCA/ADMIXTURE results are still available."
    else
        warn "visualize_results.py not found. Skipping plots."
    fi
}

print_summary() {
    log ""
    log "GAI pipeline (gnomAD v${GNOMAD_VERSION} HGDP+1kGP reference) complete."
    log "Results:"
    log "  QC report: ${RESULTS}/qc_report.txt"
    log "  PCA: ${RESULTS}/pca/merged_pca.eigenvec"
    log "  ADMIXTURE: ${RESULTS}/admixture/*.Q"
    log "  CV errors: ${RESULTS}/admixture/cv_errors.txt"
    log "  Plots: ${RESULTS}/plots/"
    log ""
}

main() {
    log "Starting GAI pipeline (gnomAD v${GNOMAD_VERSION} HGDP+1kGP reference)"
    log "Project: ${PROJECT_DIR}"
    log "Threads: ${THREADS}"

    CONDA_EXE="$(which conda || echo "")"
    if [ -n "$CONDA_EXE" ]; then
        CONDA_BASE="$(dirname "$(dirname "$CONDA_EXE")")"
        source "${CONDA_BASE}/etc/profile.d/conda.sh" || true
    fi

    step0_setup
    conda activate "${CONDA_ENV}" || true

    step1_convert
    step2_qc
    step2b_ld_prune_study
    step3_download_reference
    step4_prepare_reference
    step5_merge_prune
    step6_pca
    step7_admixture
    step8_visualize
    print_summary
}

main "$@"
