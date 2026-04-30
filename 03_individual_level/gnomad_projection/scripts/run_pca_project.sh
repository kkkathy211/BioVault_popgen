#!/usr/bin/env bash
# Path A: PCA projection of study samples onto gnomAD HGDP+1kGP reference
# space using gnomAD's pre-computed PCA loadings.
#
# Skips per-chr VCF download entirely. Total runtime ~minutes (after one-time
# Hail install of ~5-10 min). Independent of the main pipeline — safe to run
# in a separate terminal while run_pipeline_gnomad.sh is downloading chrs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"
WORKING="${SCRIPT_DIR}/working"
RESULTS="${SCRIPT_DIR}/results"
CONDA_ENV="ancestry_pipeline"

OUT_DIR="${RESULTS}/pca_projection"
mkdir -p "${OUT_DIR}"

log()  { echo "[$(date '+%H:%M:%S')] $*"; }
die()  { echo "[$(date '+%H:%M:%S')] ERROR: $*" >&2; exit 1; }

# --- Activate conda env ------------------------------------------------------
CONDA_EXE="$(which conda || echo "")"
[ -n "$CONDA_EXE" ] || die "conda not found"
CONDA_BASE="$(dirname "$(dirname "$CONDA_EXE")")"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

# --- Java check (Hail needs Java 8 or 11) -----------------------------------
if ! command -v java >/dev/null 2>&1; then
    log "Installing OpenJDK 11 into conda env (Hail dependency)..."
    conda install -n "${CONDA_ENV}" -c conda-forge openjdk=11 -y
fi

# --- Hail install (one-time, ~5-10 min) -------------------------------------
if ! python -c "import hail" 2>/dev/null; then
    log "Installing Hail (one-time, ~5-10 min)..."
    # pip is more reliable than conda for Hail; bioconda's hail can lag releases
    pip install --quiet hail
fi

# --- Inputs ------------------------------------------------------------------
STUDY_PREFIX="${WORKING}/study_qc"
for ext in bed bim fam; do
    [ -f "${STUDY_PREFIX}.${ext}" ] || die "Missing ${STUDY_PREFIX}.${ext} -- run step1+step2 of run_pipeline_gnomad.sh first"
done

log "Study PLINK: ${STUDY_PREFIX}"
log "Output dir:  ${OUT_DIR}"

# --- Run projection ----------------------------------------------------------
python "${SCRIPTS_DIR}/pca_project.py" "${STUDY_PREFIX}" "${OUT_DIR}"

# --- Overlay study on gnomAD reference PCs -----------------------------------
META_TSV="${SCRIPT_DIR}/reference/hgdp_tgp_sample_meta.tsv"
if [ -f "${META_TSV}" ]; then
    log "Building reference-overlay plot..."
    python "${SCRIPTS_DIR}/pca_overlay_plot.py" \
        "${META_TSV}" \
        "${OUT_DIR}/study_pca_projection.tsv" \
        "${OUT_DIR}/pca_overlay_pc1_pc2.png" 1 2
    python "${SCRIPTS_DIR}/pca_overlay_plot.py" \
        "${META_TSV}" \
        "${OUT_DIR}/study_pca_projection.tsv" \
        "${OUT_DIR}/pca_overlay_pc3_pc4.png" 3 4

    log "Building sub-population overlay plots..."
    python "${SCRIPTS_DIR}/pca_overlay_subpop.py" \
        "${META_TSV}" \
        "${OUT_DIR}/study_pca_projection.tsv" \
        "${OUT_DIR}/pca_subpop_all.png" 1 2
    for REGION in AFR AMR CSA EAS EUR MID OCE; do
        python "${SCRIPTS_DIR}/pca_overlay_subpop.py" \
            "${META_TSV}" \
            "${OUT_DIR}/study_pca_projection.tsv" \
            "${OUT_DIR}/pca_subpop_${REGION}.png" 1 2 "${REGION}"
    done
else
    log "WARN: ${META_TSV} not found -- skipping overlay plot."
fi

log ""
log "Done. Results:"
log "  ${OUT_DIR}/study_pca_projection.tsv"
log "  ${OUT_DIR}/pca_projection.png            (study only)"
log "  ${OUT_DIR}/pca_overlay_pc1_pc2.png       (with reference)"
log "  ${OUT_DIR}/pca_overlay_pc3_pc4.png       (with reference)"
