"""
06_AIMs_dendrogram.py
=====================
Goal: pull TT and Bermuda apart on a principled set of ancestry-informative
markers (AIMs), without cherry-picking SNPs from the Caribbean cohorts
themselves (which would be circular).

Strategy
--------
Pick AIMs entirely from gnomAD reference-population contrasts:
  Panel A: top 200 SNPs by |AF_AFR − AF_NFE|  (African vs European axis)
  Panel B: top 200 SNPs by |AF_AFR − AF_SAS|  (African vs South-Asian axis)

Then project each Caribbean cohort onto those panels and:
  1) clustered heatmap (rows = AIMs, columns = islands + the relevant gnomAD
     reference pops) with hierarchical-clustering dendrograms on both axes;
  2) PCA over the AIMs (rows = cohorts, dims = AIMs) → 2D scatter showing the
     placement of TT and Bermuda relative to AFR / NFE / SAS reference vertices.

Outputs
-------
  ../data/aims/aims_AFR_NFE.tsv         locus_key + island/ref AFs + |ΔAF|, rank
  ../data/aims/aims_AFR_SAS.tsv         "
  ../data/aims/aims_combined.tsv        union of both panels
  ../plots/aims_AFR_NFE_clustermap.pdf  clustered heatmap + dendrograms
  ../plots/aims_AFR_NFE_clustermap.png
  ../plots/aims_AFR_SAS_clustermap.pdf
  ../plots/aims_AFR_SAS_clustermap.png
  ../plots/aims_combined_pca.pdf        PCA of cohorts on combined AIMs
  ../plots/aims_combined_pca.png

Pre-filters
-----------
  * AF_AFR + AF_NFE (or AF_AFR + AF_SAS) both informative: min(AF, 1-AF) ≥ 0.05
    (drops near-monomorphic sites).
  * |ΔAF| sorted; ties broken by smaller pooled AF (slightly favours rarer,
    more diagnostic variants).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage

BASE      = Path(__file__).resolve().parents[1]
MASTER    = BASE / "data" / "master_af_table.tsv"
AIM_DIR   = BASE / "data" / "aims"
PLOTS_DIR = BASE / "plots"
AIM_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

ISLANDS     = ["BVI", "TT", "Bahamas", "Barbados", "Bermuda", "StLucia"]
GNOMAD_REFS = ["gnomAD_global", "gnomAD_AFR", "gnomAD_NFE", "gnomAD_SAS"]
TOP_N       = 200
MAF_FLOOR   = 0.05

DISPLAY = {
    "BVI": "BVI", "TT": "Trinidad &\nTobago", "Bahamas": "Bahamas",
    "Barbados": "Barbados", "Bermuda": "Bermuda", "StLucia": "St. Lucia",
    "gnomAD_global": "gnomAD global", "gnomAD_AFR": "gnomAD AFR",
    "gnomAD_NFE": "gnomAD NFE", "gnomAD_SAS": "gnomAD SAS",
}
COLORS = {
    "BVI": "#E69F00", "TT": "#56B4E9", "Bahamas": "#009E73",
    "Barbados": "#CC79A7", "Bermuda": "#D55E00", "StLucia": "#F0E442",
    "gnomAD_AFR": "#000000", "gnomAD_NFE": "#7F7F7F",
    "gnomAD_SAS": "#999999", "gnomAD_global": "#444444",
}


def maf(x: pd.Series) -> pd.Series:
    return x.clip(0, 1).combine(1 - x, min)


def pick_aims(df: pd.DataFrame, ref_a: str, ref_b: str, top_n: int) -> pd.DataFrame:
    common = (maf(df[ref_a]) >= MAF_FLOOR) & (maf(df[ref_b]) >= MAF_FLOOR)
    sub = df.loc[common].copy()
    sub["delta_abs"] = (sub[ref_a] - sub[ref_b]).abs()
    sub["pooled_af"] = (sub[ref_a] + sub[ref_b]) / 2
    sub = (sub.sort_values(["delta_abs", "pooled_af"], ascending=[False, True])
              .head(top_n)
              .reset_index(drop=True))
    sub["rank"] = sub.index + 1
    return sub


def write_panel(panel: pd.DataFrame, ref_a: str, ref_b: str, name: str) -> Path:
    cols = (["locus_key", "rsid", "rank", "delta_abs", "pooled_af"]
            + ISLANDS + GNOMAD_REFS)
    f = AIM_DIR / f"{name}.tsv"
    panel[cols].to_csv(f, sep="\t", index=False, float_format="%.6f")
    return f


# ── Heatmap with dendrograms ─────────────────────────────────────────────────
def plot_clustermap(panel: pd.DataFrame, refs_in_view: list[str],
                    label_a: str, label_b: str, fname_stem: str) -> None:
    cols = ISLANDS + refs_in_view
    M = panel.set_index("locus_key")[cols].copy()
    M.columns = [DISPLAY[c] for c in cols]

    col_colors = [COLORS[c] for c in cols]

    g = sns.clustermap(
        M,
        cmap="RdBu_r",
        center=0.5,
        vmin=0, vmax=1,
        figsize=(6.0, 9.0),
        method="average",
        metric="euclidean",
        col_cluster=True,
        row_cluster=True,
        col_colors=[col_colors],
        cbar_kws={"label": "Allele frequency"},
        xticklabels=True,
        yticklabels=False,
        dendrogram_ratio=(0.18, 0.12),
        colors_ratio=0.025,
        linewidths=0,
    )
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel(f"AIMs (top {TOP_N} {label_a}↔{label_b}, |ΔAF|-ranked)",
                            fontsize=8)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    g.fig.suptitle(
        f"AIMs from gnomAD {label_a} vs {label_b} — applied to Caribbean cohorts",
        fontsize=9, y=1.0,
    )
    pdf = PLOTS_DIR / f"{fname_stem}.pdf"
    png = PLOTS_DIR / f"{fname_stem}.png"
    g.savefig(pdf, bbox_inches="tight")
    g.savefig(png, bbox_inches="tight", dpi=300)
    plt.close(g.fig)
    print(f"[06] wrote {pdf}")
    print(f"[06] wrote {png}")


# ── PCA on combined AIMs ─────────────────────────────────────────────────────
def pca_on_aims(combined: pd.DataFrame) -> None:
    cohorts = ISLANDS + ["gnomAD_AFR", "gnomAD_NFE", "gnomAD_SAS", "gnomAD_global"]
    X = combined[cohorts].to_numpy().T            # rows = cohorts, cols = AIMs
    X = X - X.mean(axis=0, keepdims=True)
    # SVD-based PCA (avoids sklearn dependency)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    PC = U * S                                    # (n_cohorts, n_pcs)
    var_explained = (S ** 2) / np.sum(S ** 2) * 100

    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    for i, c in enumerate(cohorts):
        is_ref = c.startswith("gnomAD")
        ax.scatter(PC[i, 0], PC[i, 1],
                   s=110 if not is_ref else 80,
                   color=COLORS.get(c, "#444"),
                   edgecolor="black", linewidth=0.6,
                   marker="o" if not is_ref else "s",
                   label=DISPLAY[c].replace("\n", " "),
                   zorder=3)
        ax.annotate(DISPLAY[c].replace("\n", " "),
                    (PC[i, 0], PC[i, 1]),
                    xytext=(5, 5), textcoords="offset points",
                    fontsize=6.5, zorder=4)
    ax.axhline(0, color="grey", lw=0.4, ls="--", zorder=1)
    ax.axvline(0, color="grey", lw=0.4, ls="--", zorder=1)
    ax.set_xlabel(f"PC1 ({var_explained[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_explained[1]:.1f}%)")
    ax.set_title(
        f"PCA of cohorts on combined AIMs\n"
        f"(top {TOP_N} AFR↔NFE ∪ top {TOP_N} AFR↔SAS, deduplicated)",
        fontsize=8,
    )
    ax.tick_params(labelsize=7)
    fig.tight_layout()
    pdf = PLOTS_DIR / "aims_combined_pca.pdf"
    png = PLOTS_DIR / "aims_combined_pca.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[06] wrote {pdf}")
    print(f"[06] wrote {png}")


def main() -> None:
    master = pd.read_csv(MASTER, sep="\t")
    print(f"[06] master rows: {len(master):,}")

    panel_ne = pick_aims(master, "gnomAD_AFR", "gnomAD_NFE", TOP_N)
    panel_sa = pick_aims(master, "gnomAD_AFR", "gnomAD_SAS", TOP_N)
    write_panel(panel_ne, "gnomAD_AFR", "gnomAD_NFE", "aims_AFR_NFE")
    write_panel(panel_sa, "gnomAD_AFR", "gnomAD_SAS", "aims_AFR_SAS")
    print(f"[06] AFR/NFE panel: {len(panel_ne)}; AFR/SAS panel: {len(panel_sa)}")

    plot_clustermap(panel_ne,
                    refs_in_view=["gnomAD_AFR", "gnomAD_NFE", "gnomAD_global"],
                    label_a="AFR", label_b="NFE",
                    fname_stem="aims_AFR_NFE_clustermap")
    plot_clustermap(panel_sa,
                    refs_in_view=["gnomAD_AFR", "gnomAD_SAS", "gnomAD_global"],
                    label_a="AFR", label_b="SAS",
                    fname_stem="aims_AFR_SAS_clustermap")

    combined = pd.concat([panel_ne, panel_sa]).drop_duplicates(subset="locus_key")
    combined.to_csv(AIM_DIR / "aims_combined.tsv", sep="\t", index=False,
                    float_format="%.6f")
    print(f"[06] combined AIMs (dedup): {len(combined)}")
    pca_on_aims(combined)


if __name__ == "__main__":
    main()
