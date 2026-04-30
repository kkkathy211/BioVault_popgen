#!/usr/bin/env bash
# 02_fetch_gnomad_v4.sh   (now actually fetches gnomAD v3.1.2)
# ============================================================
# Slice gnomAD v3.1.2 main genomes sites VCFs over our pre-filtered Caribbean
# loci, in parallel by chromosome.  v3.1.2 publishes per-population AF as
# INFO fields directly:
#     AF          (global, all gnomAD genomes)
#     AF_afr      (African / African-American)
#     AF_nfe      (Non-Finnish European)
#     AF_sas      (South Asian)
# HGDP+1KGP samples are part of this cohort; the main release does not break
# them out into a separate AF_hgdp_* subset, so we use the population AFs
# directly (HGDP+1KGP is a subset of these and the population frequencies are
# very close).
#
# Public URLs:
#   https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr<N>.vcf.bgz
#
# Input
# -----
#   ../data/caribbean_loci.bed           (built by 01)
#
# Output
# ------
#   ../data/gnomad_v4_raw/chr<N>.tsv     (kept name for downstream stability)
#       columns:  CHROM POS REF ALT AF AF_afr AF_nfe AF_sas
#   ../logs/02_fetch_gnomad_v4.log
#
# Speed
# -----
#   * PARALLEL chromosomes (default 8 in flight)
#   * Pre-filtered BED (~300 k SNPs)
#   Expected wall time: 20–40 minutes on a normal connection.
#   Resume-friendly: any non-empty existing chr<N>.tsv is skipped.

set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="$(dirname "$SCRIPT_DIR")"
BED="$BASE/data/caribbean_loci.bed"
OUTDIR="$BASE/data/gnomad_v4_raw"
LOG="$BASE/logs/02_fetch_gnomad_v4.log"

mkdir -p "$OUTDIR"
: > "$LOG"

GNOMAD_BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes"
PARALLEL=${PARALLEL:-8}

QUERY_FMT='%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_afr\t%INFO/AF_nfe\t%INFO/AF_sas\n'

fetch_one() {
    local chr="$1"
    local out="$OUTDIR/chr${chr}.tsv"
    if [[ -s "$out" ]]; then
        echo "[02] chr${chr}: cached, skipping" >> "$LOG"
        return 0
    fi
    local url="$GNOMAD_BASE/gnomad.genomes.v3.1.2.sites.chr${chr}.vcf.bgz"
    local bed_chr="$OUTDIR/.chr${chr}.bed"
    awk -v c="chr${chr}" '$1==c' "$BED" > "$bed_chr"
    local n
    n=$(wc -l < "$bed_chr" | tr -d ' ')
    echo "[02] chr${chr}: $n loci → starting" >> "$LOG"
    if [[ "$n" -eq 0 ]]; then
        : > "$out"
        rm -f "$bed_chr"
        return 0
    fi
    bcftools view -R "$bed_chr" "$url" 2>>"$LOG" \
        | bcftools query -f "$QUERY_FMT" 2>>"$LOG" \
        > "$out".tmp
    if [[ -s "$out".tmp ]] || [[ "$n" -eq 0 ]]; then
        mv "$out".tmp "$out"
        local rows
        rows=$(wc -l < "$out" | tr -d ' ')
        echo "[02] chr${chr}: wrote $rows rows" >> "$LOG"
    else
        echo "[02] chr${chr}: FAILED (no output)" >> "$LOG"
        rm -f "$out".tmp
    fi
    rm -f "$bed_chr"
}

export -f fetch_one
export GNOMAD_BASE BED OUTDIR LOG QUERY_FMT

CHRS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
echo "[02] launching parallel=$PARALLEL across chrs: $CHRS" | tee -a "$LOG"

printf '%s\n' $CHRS | xargs -n1 -P "$PARALLEL" -I{} bash -c 'fetch_one "$@"' _ {}

echo "[02] all chromosomes done." | tee -a "$LOG"
ls -la "$OUTDIR" | tail -30 | tee -a "$LOG"
