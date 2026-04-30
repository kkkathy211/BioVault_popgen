#!/usr/bin/env bash
# =============================================================================
# Master runner — Caribbean FST pipeline
# Usage: PYTHON=python3.9 bash scripts/run_pipeline.sh
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON="${PYTHON:-python3}"

echo "============================================="
echo "  Caribbean Island FST Pipeline"
echo "============================================="

echo ""
echo "[Step 1] Load & merge allele frequency files …"
$PYTHON "$SCRIPT_DIR/01_load_merge.py"

echo ""
echo "[Step 2] Compute pairwise FST (Weir & Cockerham 1984) …"
$PYTHON "$SCRIPT_DIR/02_compute_fst.py"

echo ""
echo "[Step 3] Generate visualisations …"
$PYTHON "$SCRIPT_DIR/03_visualize.py"

echo ""
echo "============================================="
echo "  Pipeline complete."
echo "  Outputs:"
echo "    FST matrix : $BASE_DIR/data/fst/fst_matrix.tsv"
echo "    Plots      : $BASE_DIR/plots/"
echo "    PCA coords : $BASE_DIR/data/pca/population_pca.tsv"
echo "============================================="
