"""
Step 3 — Visualise FST results.

Outputs
-------
  plots/fst_heatmap.png           — annotated FST heatmap
  plots/fst_dendrogram.png        — hierarchical clustering tree
  plots/fst_heatmap_clustermap.png— combined clustered heatmap + dendrogram
  plots/pca_populations.png       — PCA of populations in allele-freq space
  data/pca/population_pca.tsv     — PC coordinates + variance explained

Notes
-----
* The FST matrix is symmetric with zero diagonal.  We pass 1-FST as a
  *distance* to the hierarchical clustering, so close populations
  (low FST) cluster together.
* Population-level PCA operates on the SNP × population allele frequency
  matrix (transposed for sklearn: populations are observations, SNPs are
  features).  This is a completely different thing from individual-level
  PCA — here each "individual" is an entire island, described by its
  vector of allele frequencies over ~1 M SNPs.
"""

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR  = Path(__file__).resolve().parents[1]
FST_PATH  = BASE_DIR / "data" / "fst" / "fst_matrix.tsv"
FREQ_PATH = BASE_DIR / "data" / "merged" / "merged_allele_freq.tsv"
PCA_DIR   = BASE_DIR / "data" / "pca"
PLOTS_DIR = BASE_DIR / "plots"
LOG_DIR   = BASE_DIR / "logs"

# ── Style ─────────────────────────────────────────────────────────────────────
PALETTE     = "YlOrRd"          # perceptual heatmap for FST (0 = similar)
ACCENT      = "#2c7bb6"         # colour for dendrogram + PCA
FIG_DPI     = 150

# ── Logging ───────────────────────────────────────────────────────────────────
LOG_DIR.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "03_visualize.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# 1.  FST Heatmap
# ─────────────────────────────────────────────────────────────────────────────

def plot_fst_heatmap(fst: pd.DataFrame, out: Path):
    """Plain annotated FST heatmap (no reordering)."""
    mask = np.eye(len(fst), dtype=bool)   # mask diagonal

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(
        fst,
        annot=True, fmt=".4f",
        cmap=PALETTE,
        linewidths=0.5, linecolor="white",
        mask=mask,
        ax=ax,
        cbar_kws={"label": "Pairwise FST (WC84)", "shrink": 0.8},
        vmin=0,
    )
    ax.set_title("Pairwise FST — Caribbean Island Populations\n"
                 "(Weir & Cockerham 1984)", fontsize=13, fontweight="bold", pad=12)
    ax.tick_params(axis="both", labelsize=10)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
    plt.tight_layout()
    plt.savefig(out, dpi=FIG_DPI)
    plt.close()
    log.info(f"Saved FST heatmap → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Hierarchical clustering — standalone dendrogram
# ─────────────────────────────────────────────────────────────────────────────

def plot_dendrogram(fst: pd.DataFrame, out: Path):
    """
    Use FST directly as a genetic distance matrix (not 1 − FST).
    Because all pairwise FST values here are small (0.003–0.03), using
    1 − FST as distance would compress everything to ~0.97–1.00, losing all
    visual resolution.  FST itself is already a meaningful genetic distance:
    0 = identical, higher = more differentiated.

    Run UPGMA (average-linkage) clustering — the standard in pop-gen trees.
    """
    pops = list(fst.columns)
    dist = np.clip(fst.values, 0, None)         # FST ≥ 0; diagonal = 0
    np.fill_diagonal(dist, 0.0)
    condensed = squareform(dist, checks=False)

    # UPGMA (average linkage) — most common in population genetics trees
    Z = linkage(condensed, method="average")

    fig, ax = plt.subplots(figsize=(8, 5))
    dendrogram(
        Z,
        labels=pops,
        ax=ax,
        color_threshold=0.6 * max(Z[:, 2]),
        above_threshold_color="grey",
        leaf_font_size=12,
        link_color_func=lambda k: ACCENT,
    )
    ax.set_title("UPGMA Clustering of Caribbean Populations\n"
                 "Distance = FST (WC84 pairwise)", fontsize=13, fontweight="bold")
    ax.set_xlabel("Population", fontsize=11)
    ax.set_ylabel("FST distance", fontsize=11)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.4f"))
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    plt.savefig(out, dpi=FIG_DPI)
    plt.close()
    log.info(f"Saved dendrogram → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Combined clustered heatmap + dendrogram  (seaborn clustermap)
# ─────────────────────────────────────────────────────────────────────────────

def plot_clustermap(fst: pd.DataFrame, out: Path):
    """
    Seaborn clustermap: reorders rows/columns by hierarchical clustering and
    draws marginal dendrograms.  We pre-compute the linkage from condensed FST
    distances so seaborn doesn't try to re-interpret the symmetric matrix.
    """
    dist_condensed = squareform(np.clip(fst.values, 0, None), checks=False)
    Z = linkage(dist_condensed, method="average")

    g = sns.clustermap(
        fst,
        row_linkage=Z,
        col_linkage=Z,
        cmap=PALETTE,
        annot=True, fmt=".4f",
        linewidths=0.5, linecolor="white",
        figsize=(8, 7),
        cbar_pos=(0.02, 0.8, 0.03, 0.18),
        dendrogram_ratio=0.2,
        vmin=0,
    )
    g.fig.suptitle("Clustered FST Heatmap — Caribbean Islands",
                   y=1.02, fontsize=13, fontweight="bold")
    g.ax_heatmap.set_xticklabels(
        g.ax_heatmap.get_xticklabels(), rotation=30, ha="right", fontsize=10)
    g.ax_heatmap.set_yticklabels(
        g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=10)
    plt.savefig(out, dpi=FIG_DPI, bbox_inches="tight")
    plt.close()
    log.info(f"Saved clustermap → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Population-level PCA on allele frequency matrix
# ─────────────────────────────────────────────────────────────────────────────

def plot_population_pca(freq_path: Path, pca_out_dir: Path, plot_out: Path):
    """
    PCA on the allele-frequency matrix.

    * Rows   = SNPs  (~1 M)
    * Columns = populations  (6)
    → Transpose so PCA sees 6 "observations" (populations) × 1 M "features" (SNPs).

    After mean-imputation of any residual NaNs and optional standardisation,
    we compute all possible PCs (min(6,1M) = 6, but only 5 are non-trivial
    because the data matrix has 6 samples → at most 5 informative PCs).

    This answers: "which populations have similar genome-wide allele freq profiles?"
    It is *not* individual-level ancestry.
    """
    log.info("Running population-level PCA on allele frequencies …")
    freq = pd.read_csv(freq_path, sep="\t", index_col="locus_key")

    pops  = list(freq.columns)
    n_pops = len(pops)

    # X: populations × SNPs   (each population is a row, each SNP a feature)
    X = freq.values.T.astype(float)          # shape: (6, ~1M)

    # Impute any residual NaN with column (SNP) mean
    imp = SimpleImputer(strategy="mean")
    X   = imp.fit_transform(X.T).T           # impute along SNP axis

    # Standardise each SNP (subtract mean, divide by std) — optional but
    # recommended so rare-variant SNPs don't dominate.
    scaler = StandardScaler()
    X_sc   = scaler.fit_transform(X.T).T     # still (6, S)

    n_components = min(n_pops - 1, n_pops)   # max meaningful PCs for n samples
    pca  = PCA(n_components=n_components)
    scores = pca.fit_transform(X_sc)         # (6, n_components)

    var_exp = pca.explained_variance_ratio_ * 100

    # Save PC coordinates
    pca_out_dir.mkdir(parents=True, exist_ok=True)
    pc_df = pd.DataFrame(scores, index=pops,
                          columns=[f"PC{i+1}" for i in range(n_components)])
    pc_df.index.name = "population"
    pc_df.to_csv(pca_out_dir / "population_pca.tsv", sep="\t", float_format="%.6f")

    eigenval_path = pca_out_dir / "population_pca_eigenval.tsv"
    pd.Series(pca.explained_variance_ratio_, name="var_explained").to_csv(
        eigenval_path, sep="\t")

    # ── Plot ──────────────────────────────────────────────────────────────
    colours = plt.cm.get_cmap("tab10", n_pops)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for ax_idx, (pc_x, pc_y) in enumerate([("PC1", "PC2"), ("PC2", "PC3")]):
        ax = axes[ax_idx]
        if pc_y not in pc_df.columns:
            ax.set_visible(False)
            continue

        for k, pop in enumerate(pops):
            ax.scatter(pc_df.loc[pop, pc_x], pc_df.loc[pop, pc_y],
                       color=colours(k), s=120, edgecolors="k",
                       linewidths=0.7, zorder=3, label=pop)
            ax.annotate(pop, (pc_df.loc[pop, pc_x], pc_df.loc[pop, pc_y]),
                        textcoords="offset points", xytext=(7, 4),
                        fontsize=9, color="dimgray")

        pcx_idx = int(pc_x[2:]) - 1
        pcy_idx = int(pc_y[2:]) - 1
        ax.set_xlabel(f"{pc_x} ({var_exp[pcx_idx]:.1f}% var)", fontsize=11)
        ax.set_ylabel(f"{pc_y} ({var_exp[pcy_idx]:.1f}% var)", fontsize=11)
        ax.axhline(0, color="lightgray", lw=0.8, ls="--")
        ax.axvline(0, color="lightgray", lw=0.8, ls="--")
        ax.grid(True, alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(fontsize=8, loc="best", framealpha=0.7)

    fig.suptitle("Population-Level PCA — Allele Frequency Space\n"
                 "(each point = one island population)",
                 fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(plot_out, dpi=FIG_DPI)
    plt.close()
    log.info(f"Saved population PCA → {plot_out}")

    for i, v in enumerate(var_exp, 1):
        log.info(f"  PC{i}: {v:.2f}% variance explained")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    log.info("Loading FST matrix …")
    fst = pd.read_csv(FST_PATH, sep="\t", index_col=0)
    log.info(f"  Populations: {list(fst.columns)}")

    plot_fst_heatmap(fst, PLOTS_DIR / "fst_heatmap.png")
    plot_dendrogram (fst, PLOTS_DIR / "fst_dendrogram.png")
    plot_clustermap (fst, PLOTS_DIR / "fst_heatmap_clustermap.png")
    plot_population_pca(FREQ_PATH, PCA_DIR, PLOTS_DIR / "pca_populations.png")

    log.info("Step 3 complete.  All plots saved to plots/")


if __name__ == "__main__":
    main()
