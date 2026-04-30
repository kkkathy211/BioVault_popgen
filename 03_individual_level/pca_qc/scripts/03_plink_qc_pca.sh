#!/usr/bin/env bash
# =============================================================================
# Step 3: PLINK QC → LD pruning → PCA
#
# Requirements: PLINK 1.9  (brew install plink  OR  conda install -c bioconda plink)
#
# Usage: bash 03_plink_qc_pca.sh
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
PLINK_DIR="$BASE_DIR/data/plink"
PCA_DIR="$BASE_DIR/data/pca"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$PCA_DIR" "$LOG_DIR"

INPUT="$PLINK_DIR/genotypes"          # .ped + .map
QC1="$PLINK_DIR/qc_pass"
PRUNED="$PLINK_DIR/pruned"
PCA_OUT="$PCA_DIR/pca"

# ---------------------------------------------------------------------------
# QC thresholds (BioVault paper values)
# ---------------------------------------------------------------------------
GENO=0.05       # max per-SNP missingness   (call rate > 95 %)
MIND=0.1        # max per-individual missingness
MAF=0.01        # minor allele frequency
HWE=1e-4        # Hardy-Weinberg equilibrium p-value threshold

# LD pruning parameters
LD_WINDOW=50
LD_STEP=5
LD_R2=0.2

# ---------------------------------------------------------------------------
# 1. Convert PED/MAP → binary BED/BIM/FAM, apply QC filters
# ---------------------------------------------------------------------------
echo "[$(date)] Step 3a: QC filtering …"
plink \
    --file     "$INPUT" \
    --geno     "$GENO" \
    --mind     "$MIND" \
    --maf      "$MAF" \
    --hwe      "$HWE" \
    --make-bed \
    --out      "$QC1" \
    --allow-no-sex \
    2>&1 | tee "$LOG_DIR/plink_qc.log"

echo "[$(date)] QC complete. Output: $QC1.{bed,bim,fam}"

# ---------------------------------------------------------------------------
# 2. LD pruning — generate list of pruned-in SNPs
# ---------------------------------------------------------------------------
echo "[$(date)] Step 3b: LD pruning …"
plink \
    --bfile    "$QC1" \
    --indep-pairwise "$LD_WINDOW" "$LD_STEP" "$LD_R2" \
    --out      "$PRUNED" \
    --allow-no-sex \
    2>&1 | tee "$LOG_DIR/plink_prune.log"

echo "[$(date)] LD pruning complete. Pruned-in list: $PRUNED.prune.in"

# ---------------------------------------------------------------------------
# 3. Extract pruned SNPs and run PCA (default: 20 PCs)
# ---------------------------------------------------------------------------
echo "[$(date)] Step 3c: Running PCA …"
plink \
    --bfile    "$QC1" \
    --extract  "$PRUNED.prune.in" \
    --pca      20 \
    --out      "$PCA_OUT" \
    --allow-no-sex \
    2>&1 | tee "$LOG_DIR/plink_pca.log"

echo "[$(date)] PCA complete."
echo "  Eigenvectors : $PCA_OUT.eigenvec"
echo "  Eigenvalues  : $PCA_OUT.eigenval"
echo ""
echo "Run: python scripts/04_plot_pca.py"
