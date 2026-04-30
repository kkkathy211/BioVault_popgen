"""
Step 3b (Python alternative): Run PCA directly in Python when PLINK is unavailable.

Applies the same QC thresholds as the PLINK shell script:
    call rate > 95 %  (GENO 0.05)
    MAF > 1 %
    HWE p > 1e-4
    LD pruning via variance inflation (approximate)

Input:  data/merged/genotype_matrix_numeric.tsv
Output: data/pca/pca.eigenvec   (same format as PLINK)
        data/pca/pca.eigenval
"""

import logging
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR  = Path(__file__).resolve().parents[1]
MERGED    = BASE_DIR / "data" / "merged"
PCA_DIR   = BASE_DIR / "data" / "pca"
LOG_DIR   = BASE_DIR / "logs"

N_PCS     = 20
GENO      = 0.05      # max per-SNP missing rate
MIND      = 0.10      # max per-individual missing rate
MAF       = 0.01
HWE_P     = 1e-4
LD_R2     = 0.2       # approximate LD pruning threshold

# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "03b_python_pca.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# QC helpers
# ---------------------------------------------------------------------------
def filter_call_rate(mat: pd.DataFrame, geno: float, mind: float) -> pd.DataFrame:
    """Remove low-call-rate SNPs (rows) and individuals (columns)."""
    snp_missing  = mat.isna().mean(axis=1)
    ind_missing  = mat.isna().mean(axis=0)
    log.info(f"  Before call rate filter: {mat.shape}")
    mat = mat.loc[snp_missing <= geno, ind_missing <= mind]
    log.info(f"  After  call rate filter: {mat.shape}")
    return mat


def filter_maf(mat: pd.DataFrame, maf_thresh: float) -> pd.DataFrame:
    """Remove SNPs with MAF below threshold."""
    # dosage coded 0/1/2; mean dosage ≈ 2 * alt freq
    alt_freq = mat.mean(axis=1, skipna=True) / 2.0
    maf_vals = alt_freq.clip(upper=0.5)           # flip if > 0.5
    maf_vals = maf_vals.where(alt_freq <= 0.5, 1 - alt_freq)
    keep = maf_vals >= maf_thresh
    log.info(f"  SNPs after MAF>{maf_thresh}: {keep.sum()} / {len(keep)}")
    return mat.loc[keep]


def filter_hwe(mat: pd.DataFrame, p_thresh: float) -> pd.DataFrame:
    """Remove SNPs failing HWE (chi-squared test on genotype counts)."""
    keep = []
    for rsid, row in mat.iterrows():
        dosage = row.dropna().values
        n_hom_ref = (dosage == 0).sum()
        n_het     = (dosage == 1).sum()
        n_hom_alt = (dosage == 2).sum()
        n = n_hom_ref + n_het + n_hom_alt
        if n == 0:
            keep.append(False)
            continue
        p_alt = (2 * n_hom_alt + n_het) / (2 * n)
        p_ref = 1 - p_alt
        exp_hom_ref = n * p_ref ** 2
        exp_het     = n * 2 * p_ref * p_alt
        exp_hom_alt = n * p_alt ** 2
        obs = [n_hom_ref, n_het, n_hom_alt]
        exp = [exp_hom_ref, exp_het, exp_hom_alt]
        # Avoid chi2 with zero expected cells
        if any(e < 1e-6 for e in exp):
            keep.append(True)
            continue
        chi2 = sum((o - e) ** 2 / e for o, e in zip(obs, exp))
        p_val = stats.chi2.sf(chi2, df=1)
        keep.append(p_val >= p_thresh)
    result = mat.loc[keep]
    log.info(f"  SNPs after HWE p>{p_thresh}: {result.shape[0]} / {mat.shape[0]}")
    return result


def ld_prune(mat: pd.DataFrame, r2_thresh: float) -> pd.DataFrame:
    """
    Greedy LD pruning: iterate SNPs; remove any that are r^2 > threshold
    with an already-selected SNP within the same chromosome window.
    This is an approximation of PLINK's --indep-pairwise.
    """
    log.info(f"  LD pruning (r^2 threshold={r2_thresh}) …")
    imputer = SimpleImputer(strategy="mean")
    X = imputer.fit_transform(mat.values.T)   # shape: samples × SNPs
    X = X.T                                    # shape: SNPs × samples

    # Z-score standardise each SNP
    std = X.std(axis=1, keepdims=True)
    std[std == 0] = 1
    Xz = (X - X.mean(axis=1, keepdims=True)) / std

    n_snps = Xz.shape[0]
    selected = np.ones(n_snps, dtype=bool)

    for i in range(n_snps):
        if not selected[i]:
            continue
        for j in range(i + 1, min(i + 51, n_snps)):   # window of 50
            if not selected[j]:
                continue
            r = np.corrcoef(Xz[i], Xz[j])[0, 1]
            if r ** 2 > r2_thresh:
                selected[j] = False

    pruned = mat.iloc[selected]
    log.info(f"  After LD pruning: {pruned.shape[0]} SNPs retained from {n_snps}")
    return pruned


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    PCA_DIR.mkdir(parents=True, exist_ok=True)

    log.info("Loading numeric genotype matrix …")
    mat = pd.read_csv(MERGED / "genotype_matrix_numeric.tsv", sep="\t", index_col="rsid")
    log.info(f"  Loaded: {mat.shape[0]} SNPs × {mat.shape[1]} samples")

    log.info("Applying QC filters …")
    mat = filter_call_rate(mat, GENO, MIND)
    mat = filter_maf(mat, MAF)
    mat = filter_hwe(mat, HWE_P)
    mat = ld_prune(mat, LD_R2)

    # Impute remaining missing values with mean dosage (per SNP)
    log.info("Imputing residual missing values …")
    imputer = SimpleImputer(strategy="mean")
    X = imputer.fit_transform(mat.values.T)   # shape: samples × SNPs

    # PCA
    log.info(f"Running PCA (n_components={N_PCS}) …")
    n_pcs = min(N_PCS, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n_pcs)
    scores = pca.fit_transform(X)             # shape: samples × PCs

    samples = list(mat.columns)

    # Write eigenvec (PLINK format: FID IID PC1 PC2 … PCn)
    eigenvec_path = PCA_DIR / "pca.eigenvec"
    with open(eigenvec_path, "w") as f:
        for i, sid in enumerate(samples):
            pc_vals = " ".join(f"{v:.6f}" for v in scores[i])
            f.write(f"{sid} {sid} {pc_vals}\n")
    log.info(f"Saved eigenvectors → {eigenvec_path}")

    # Write eigenval
    eigenval_path = PCA_DIR / "pca.eigenval"
    with open(eigenval_path, "w") as f:
        for val in pca.explained_variance_:
            f.write(f"{val:.6f}\n")
    log.info(f"Saved eigenvalues → {eigenval_path}")

    # Variance explained
    var_exp = pca.explained_variance_ratio_ * 100
    for i, v in enumerate(var_exp[:5], 1):
        log.info(f"  PC{i}: {v:.2f}% variance explained")

    log.info("Step 3b complete. Run: python scripts/04_plot_pca.py")


if __name__ == "__main__":
    main()
