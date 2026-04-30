"""
Step 4: Plot PCA results — PC1 vs PC2 (and PC3 vs PC4 if available).

Input:  data/pca/pca.eigenvec   (PLINK format or Python output from step 3b)
        data/pca/pca.eigenval   (optional — used for axis labels)
        data/pca/sample_metadata.tsv  (optional — columns: sample_id, population/label)

Output: plots/pca_pc1_pc2.png
        plots/pca_pc3_pc4.png  (if ≥ 4 PCs available)
"""

import logging
from typing import Optional
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR   = Path(__file__).resolve().parents[1]
PCA_DIR    = BASE_DIR / "data" / "pca"
PLOTS_DIR  = BASE_DIR / "plots"
LOG_DIR    = BASE_DIR / "logs"

EIGENVEC   = PCA_DIR / "pca.eigenvec"
EIGENVAL   = PCA_DIR / "pca.eigenval"
METADATA   = PCA_DIR / "sample_metadata.tsv"   # optional

# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "04_plot.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_eigenvec(path: Path) -> pd.DataFrame:
    """Load PLINK eigenvec file → DataFrame with columns [sample_id, PC1, PC2, …]."""
    df = pd.read_csv(path, sep=r"\s+", header=None)
    n_pcs = df.shape[1] - 2
    df.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, n_pcs + 1)]
    df["sample_id"] = df["IID"].astype(str)
    return df


def load_eigenval(path: Path) -> list[float]:
    """Load eigenvalues, compute % variance explained."""
    vals = pd.read_csv(path, header=None)[0].tolist()
    total = sum(vals)
    return [v / total * 100 for v in vals] if total > 0 else vals


def assign_colors(df: pd.DataFrame, label_col: Optional[str]) -> tuple:
    """Return (colors_per_row, colormap_legend_handles)."""
    if label_col and label_col in df.columns:
        labels = df[label_col].astype(str)
        unique  = sorted(labels.unique())
        palette = cm.get_cmap("tab10", len(unique))
        color_map = {lbl: palette(i) for i, lbl in enumerate(unique)}
        colors = [color_map[l] for l in labels]
        handles = [
            plt.Line2D([0], [0], marker="o", color="w",
                       markerfacecolor=color_map[l], markersize=9, label=l)
            for l in unique
        ]
        return colors, handles
    else:
        # Colour by sample index
        n = len(df)
        palette = cm.get_cmap("tab10", n)
        colors = [palette(i) for i in range(n)]
        handles = [
            plt.Line2D([0], [0], marker="o", color="w",
                       markerfacecolor=colors[i], markersize=9,
                       label=df["sample_id"].iloc[i])
            for i in range(n)
        ]
        return colors, handles


def scatter_pca(df, pc_x, pc_y, var_exp, colors, handles, out_path):
    """Draw and save a single PCA scatter plot."""
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(
        df[pc_x], df[pc_y],
        c=colors, s=80, edgecolors="k", linewidths=0.5, alpha=0.9,
    )

    # Label each point with sample_id
    for _, row in df.iterrows():
        ax.annotate(
            row["sample_id"],
            (row[pc_x], row[pc_y]),
            textcoords="offset points", xytext=(6, 4),
            fontsize=7, color="dimgray",
        )

    pc_x_num = int(pc_x.replace("PC", ""))
    pc_y_num = int(pc_y.replace("PC", ""))

    x_label = f"{pc_x}"
    y_label = f"{pc_y}"
    if var_exp:
        x_label += f" ({var_exp[pc_x_num-1]:.1f}% var)"
        y_label += f" ({var_exp[pc_y_num-1]:.1f}% var)"

    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_title("Ancestry PCA — Mock Synthetic Data", fontsize=13, fontweight="bold")
    ax.legend(handles=handles, title="Sample", fontsize=8,
              loc="best", framealpha=0.7)
    ax.axhline(0, color="lightgray", linewidth=0.8, linestyle="--")
    ax.axvline(0, color="lightgray", linewidth=0.8, linestyle="--")
    ax.grid(True, alpha=0.25)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    log.info(f"Saved plot → {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    log.info("Loading eigenvectors …")
    df = load_eigenvec(EIGENVEC)
    log.info(f"  {len(df)} samples, {sum(1 for c in df.columns if c.startswith('PC'))} PCs")

    var_exp = []
    if EIGENVAL.exists():
        var_exp = load_eigenval(EIGENVAL)
        log.info(f"  PC1: {var_exp[0]:.1f}%  PC2: {var_exp[1]:.1f}%")

    # Optionally join population labels
    label_col = None
    if METADATA.exists():
        meta = pd.read_csv(METADATA, sep="\t", dtype={"sample_id": str})
        df["sample_id"] = df["sample_id"].astype(str)
        df = df.merge(meta, on="sample_id", how="left")
        # First non-ID column in metadata is used as the colour label
        label_col = [c for c in meta.columns if c != "sample_id"][0]
        log.info(f"  Loaded metadata; colouring by '{label_col}'")

    colors, handles = assign_colors(df, label_col)

    # PC1 vs PC2
    scatter_pca(df, "PC1", "PC2", var_exp, colors, handles,
                PLOTS_DIR / "pca_pc1_pc2.png")

    # PC3 vs PC4 (if available)
    pc_cols = [c for c in df.columns if c.startswith("PC")]
    if len(pc_cols) >= 4:
        scatter_pca(df, "PC3", "PC4", var_exp, colors, handles,
                    PLOTS_DIR / "pca_pc3_pc4.png")

    log.info("Step 4 complete.")


if __name__ == "__main__":
    main()
