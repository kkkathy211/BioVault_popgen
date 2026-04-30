"""
Step 2 — Compute pairwise FST between all island populations.

Estimator: Weir & Cockerham 1984 (WC84) ratio-of-averages
──────────────────────────────────────────────────────────
WHY FST FROM ALLELE FREQUENCIES?
    FST measures the proportion of total genetic variance that resides
    *between* populations.  For a biallelic SNP:

        FST  =  (H_T − H_S) / H_T

    where H_T = heterozygosity in the pooled "meta-population" and
    H_S = mean within-population heterozygosity.  Both quantities are
    computable directly from allele frequencies:

        H_i  = 2 · p_i · (1 − p_i)         (expected het in pop i)
        H_S  = mean(H_i)
        p̄    = weighted mean allele frequency
        H_T  = 2 · p̄ · (1 − p̄)

    The WC84 pairwise θ̂ goes further by correcting for finite sample
    size (allele_number 2n), giving an unbiased estimate even with small
    or unequal cohorts — important here because n ranges from 152 to 976.

PAIRWISE WC84 FORMULA (r = 2 populations)
    For each SNP s between populations i and j:

        n̄   = (n_i + n_j) / 2
        p̄   = (n_i · p_i + n_j · p_j) / (n_i + n_j)   [weighted mean]
        nc  = n̄ − (n_i² + n_j²) / (2 · n̄)              [size correction]
        s²  = [n_i·(p_i−p̄)² + n_j·(p_j−p̄)²] / n̄      [between-pop var]

        h_i = 2·n_i·p_i·(1−p_i) / (2·n_i − 1)          [within-pop het]
        h̄   = (h_i + h_j) / 2

        a_s = (n̄/nc) · [s² − (p̄·(1−p̄) − s²/2 − h̄/4) / (n̄−1)]
        b_s = (n̄/(n̄−1)) · [p̄·(1−p̄) − s²/2 − (2·n̄−1)·h̄/(4·n̄)]
        c_s = h̄ / 2

        θ̂_s = a_s / (a_s + b_s + c_s)          (per-SNP estimate)

    Genome-wide pairwise FST (ratio-of-averages, more stable):
        FST  = Σ a_s / Σ(a_s + b_s + c_s)

LIMITATIONS OF AGGREGATE DATA
    • Cannot correct for relatedness within a population.
    • Cannot perform Hardy–Weinberg tests (need genotype counts, not just freq).
    • Sampling error in allele_freq is unknown without individual-level data.
    • Population structure *within* an island is invisible.
    • Phase information and LD patterns are lost.

Input
-----
    data/merged/merged_allele_freq.tsv
    data/merged/merged_allele_number.tsv

Output
------
    data/fst/fst_matrix.tsv           — population × population FST
    data/fst/fst_per_snp/             — per-SNP FST for each pair (optional)
"""

import logging
import itertools
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
MERGED   = BASE_DIR / "data" / "merged"
FST_DIR  = BASE_DIR / "data" / "fst"
LOG_DIR  = BASE_DIR / "logs"

SAVE_PER_SNP = False   # set True to write per-SNP FST TSVs (large files)

# ── Logging ───────────────────────────────────────────────────────────────────
LOG_DIR.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "02_compute_fst.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Core FST functions
# ─────────────────────────────────────────────────────────────────────────────

def wc84_pairwise_components(
    p1: np.ndarray, n1: np.ndarray,
    p2: np.ndarray, n2: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute per-SNP WC84 a, b, c components for a pair of populations.
    All inputs are 1-D numpy arrays of length S (number of SNPs).

    Returns (a, b, c) arrays of shape (S,).

    Notation follows Weir & Cockerham 1984, equations 2–5:
        r       = number of populations  (2 here)
        n_i     = allele count in population i  (= 2 × individuals)
        n_total = n1 + n2
        n_bar   = n_total / r               (arithmetic mean)
        nc      = (n_total - Σ n_i²/n_total) / (r-1)   [size-correction]
        p_bar   = weighted mean allele frequency
        s²      = between-pop variance of allele freq
        h_i     = within-pop expected heterozygosity (sample-size corrected)
    """
    r       = 2.0
    n_total = n1 + n2
    n_bar   = n_total / r                              # = (n1+n2)/2

    # ── Weighted mean allele frequency ────────────────────────────────────
    p_bar = (n1 * p1 + n2 * p2) / n_total

    # ── Sample-size correction factor  nc  (WC84 eq. 2) ──────────────────
    # nc = (n_total - Σ(n_i²)/n_total) / (r-1)
    # For r=2: nc = n_total - (n1² + n2²)/n_total
    nc = (n_total - (n1**2 + n2**2) / n_total) / (r - 1)

    # ── Between-population variance  s²  (WC84 eq. 3) ────────────────────
    s2 = (n1 * (p1 - p_bar)**2 + n2 * (p2 - p_bar)**2) / ((r - 1) * n_bar)

    # ── Within-population heterozygosity, sample-size corrected ──────────
    # h_i = 2·n_i·p_i·(1-p_i) / (2·n_i - 1)
    h1 = np.where(n1 > 1, 2 * n1 * p1 * (1 - p1) / (2 * n1 - 1), np.nan)
    h2 = np.where(n2 > 1, 2 * n2 * p2 * (1 - p2) / (2 * n2 - 1), np.nan)
    h_bar = (h1 + h2) / r                             # mean across populations

    # ── WC84 numerator and denominator components ─────────────────────────
    # a  (between-population component)
    inner = (p_bar * (1 - p_bar)
             - ((r - 1) / r) * s2
             - h_bar / 4.0)
    a = (n_bar / nc) * (s2 - inner / (n_bar - 1))

    # b  (within-population component)
    b = (n_bar / (n_bar - 1)) * (
        p_bar * (1 - p_bar)
        - ((r - 1) / r) * s2
        - (2 * n_bar - 1) * h_bar / (4 * n_bar)
    )

    # c  (within-individual component — heterozygosity)
    c = h_bar / 2.0

    return a, b, c


def genome_wide_fst(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """
    Ratio-of-averages genome-wide FST from WC84 components.
    Ignores SNPs where any component is NaN.
    Returns NaN if denominator sums to zero.
    """
    valid = np.isfinite(a) & np.isfinite(b) & np.isfinite(c)
    a_sum = a[valid].sum()
    denom = (a[valid] + b[valid] + c[valid]).sum()
    return float(a_sum / denom) if denom != 0 else float("nan")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    FST_DIR.mkdir(parents=True, exist_ok=True)

    log.info("Loading merged matrices …")
    freq = pd.read_csv(MERGED / "merged_allele_freq.tsv",   sep="\t", index_col="locus_key")
    n_df = pd.read_csv(MERGED / "merged_allele_number.tsv", sep="\t", index_col="locus_key")

    populations = list(freq.columns)
    n_pops      = len(populations)
    n_snps      = len(freq)
    log.info(f"  {n_pops} populations × {n_snps:,} SNPs")

    # Pre-convert to numpy for speed
    freq_np = freq.values.astype(float)   # shape: (S, P)
    n_np    = n_df.values.astype(float)   # shape: (S, P)

    # ── Compute pairwise FST ──────────────────────────────────────────────
    fst_matrix = pd.DataFrame(np.zeros((n_pops, n_pops)),
                               index=populations, columns=populations)

    pairs = list(itertools.combinations(range(n_pops), 2))
    log.info(f"Computing {len(pairs)} pairwise FST values …")

    if SAVE_PER_SNP:
        per_snp_dir = FST_DIR / "fst_per_snp"
        per_snp_dir.mkdir(exist_ok=True)

    for i, j in pairs:
        pop_i, pop_j = populations[i], populations[j]

        p1, n1 = freq_np[:, i], n_np[:, i]
        p2, n2 = freq_np[:, j], n_np[:, j]

        a, b, c = wc84_pairwise_components(p1, n1, p2, n2)
        fst_val = genome_wide_fst(a, b, c)

        fst_matrix.loc[pop_i, pop_j] = fst_val
        fst_matrix.loc[pop_j, pop_i] = fst_val
        log.info(f"  FST({pop_i} vs {pop_j}) = {fst_val:.6f}")

        if SAVE_PER_SNP:
            per_snp = pd.DataFrame({
                "rsid"  : freq.index,
                "a"     : a, "b": b, "c": c,
                "fst_snp": np.where(a + b + c > 0, a / (a + b + c), np.nan),
            })
            per_snp.to_csv(per_snp_dir / f"fst_{pop_i}_vs_{pop_j}.tsv",
                           sep="\t", index=False)

    # Diagonal is zero by definition
    np.fill_diagonal(fst_matrix.values, 0.0)

    # ── Save matrix ───────────────────────────────────────────────────────
    out_path = FST_DIR / "fst_matrix.tsv"
    fst_matrix.to_csv(out_path, sep="\t", float_format="%.6f")
    log.info(f"\nPairwise FST matrix saved → {out_path}")
    log.info("\n" + fst_matrix.round(5).to_string())

    log.info("\nStep 2 complete.")


if __name__ == "__main__":
    main()
