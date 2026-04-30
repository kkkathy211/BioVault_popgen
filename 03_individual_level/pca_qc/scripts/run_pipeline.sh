#!/usr/bin/env bash
# =============================================================================
# Master runner — executes all pipeline steps in order.
#
# Usage:
#   bash run_pipeline.sh              # uses Python PCA (no PLINK required)
#   bash run_pipeline.sh --plink      # uses PLINK for QC + PCA (Step 3)
# =============================================================================
set -euo pipefail

USE_PLINK=false
for arg in "$@"; do
  [[ "$arg" == "--plink" ]] && USE_PLINK=true
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

echo "============================================="
echo "  Ancestry PCA Pipeline"
echo "============================================="

echo ""
PYTHON="${PYTHON:-python3}"

echo "[Step 1] Merging genotype files …"
$PYTHON "$SCRIPT_DIR/01_merge_genotypes.py"

echo ""
echo "[Step 2] Encoding genotypes → numeric + PLINK PED/MAP …"
$PYTHON "$SCRIPT_DIR/02_encode_genotypes.py"

echo ""
if $USE_PLINK; then
  echo "[Step 3] PLINK QC + PCA …"
  bash "$SCRIPT_DIR/03_plink_qc_pca.sh"
else
  echo "[Step 3b] Python QC + PCA (use --plink to switch to PLINK) …"
  $PYTHON "$SCRIPT_DIR/03b_python_pca.py"
fi

echo ""
echo "[Step 4] Plotting PCA …"
$PYTHON "$SCRIPT_DIR/04_plot_pca.py"

echo ""
echo "============================================="
echo "  Pipeline complete."
echo "  PCA plot: $BASE_DIR/plots/pca_pc1_pc2.png"
echo "============================================="
