"""
figure3_publishable.py
======================
Publication-quality Caribbean population structure figure (Nature-style).

Panels
------
a  Population PCA (PC1 × PC2, allele-frequency space)
b  Pairwise FST heat map (UPGMA-ordered, annotated)
c  UPGMA phylogenetic tree (FST as distance)

Outputs
-------
  plots/figure3_publishable.pdf  — vector PDF, fonts embedded (submission-ready)
  plots/figure3_publishable.png  — 300 dpi raster

Usage
-----
  python scripts/figure3_publishable.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE          = Path(__file__).resolve().parents[1]
FST_PATH      = BASE / "data" / "fst"    / "fst_matrix.tsv"
PCA_PATH      = BASE / "data" / "pca"    / "population_pca.tsv"
EIGENVAL_PATH = BASE / "data" / "pca"    / "population_pca_eigenval.tsv"
PLOTS_DIR     = BASE / "plots"
OUT_PDF       = PLOTS_DIR / "figure3_publishable.pdf"
OUT_PNG       = PLOTS_DIR / "figure3_publishable.png"

# ── Publication-grade rcParams ─────────────────────────────────────────────────
plt.rcParams.update({
    # Fonts
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
    "font.size":          7,
    "axes.labelsize":     7,
    "axes.titlesize":     7.5,
    "xtick.labelsize":    6.5,
    "ytick.labelsize":    6.5,
    "legend.fontsize":    6.5,
    # Lines / ticks
    "axes.linewidth":     0.6,
    "xtick.major.width":  0.6,
    "ytick.major.width":  0.6,
    "xtick.major.size":   2.5,
    "ytick.major.size":   2.5,
    # PDF export: embed fonts (required by most journals)
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
})

# ── Population metadata ────────────────────────────────────────────────────────
# Short keys → display names (two-line where needed to save x-axis space)
DISPLAY = {
    "BVI":     "BVI",
    "TT":      "Trinidad &\nTobago",
    "Bahamas": "Bahamas",
    "Barbados": "Barbados",
    "Bermuda": "Bermuda",
    "StLucia": "St. Lucia",
}

# Okabe-Ito colorblind-safe palette
COLORS = {
    "BVI":     "#E69F00",   # orange
    "TT":      "#56B4E9",   # sky blue
    "Bahamas": "#009E73",   # green
    "Barbados": "#CC79A7",   # pink/mauve
    "Bermuda": "#0072B2",   # dark blue
    "StLucia": "#D55E00",   # vermilion
}

# ── Load data ──────────────────────────────────────────────────────────────────
fst      = pd.read_csv(FST_PATH,      sep="\t", index_col=0)
pca      = pd.read_csv(PCA_PATH,      sep="\t", index_col="population")
eigenval = pd.read_csv(EIGENVAL_PATH, sep="\t", index_col=0)

pops    = list(fst.columns)                              # original order
var_exp = eigenval["var_explained"].values * 100         # % variance per PC

# ── Hierarchical clustering (UPGMA on FST) ────────────────────────────────────
dist = np.clip(fst.values, 0.0, None)
np.fill_diagonal(dist, 0.0)
condensed = squareform(dist, checks=False)
Z         = linkage(condensed, method="average")         # UPGMA
ord_pops  = [pops[i] for i in leaves_list(Z)]           # UPGMA leaf order

# ── Figure layout ──────────────────────────────────────────────────────────────
# ~7.4 × 3.4 in  ≈  188 × 86 mm  (Nature two-column figure width)
FIG_W, FIG_H = 7.4, 3.4
fig = plt.figure(figsize=(FIG_W, FIG_H))

gs = gridspec.GridSpec(
    1, 3,
    figure=fig,
    width_ratios=[1.10, 1.10, 0.85],
    wspace=0.46,
    left=0.07, right=0.97,
    top=0.88, bottom=0.16,
)

def _panel_label(ax, letter):
    """Bold upper-left panel label (Nature convention)."""
    ax.text(
        -0.16, 1.12, letter,
        transform=ax.transAxes,
        fontsize=10, fontweight="bold",
        va="top", ha="left",
    )


# ══════════════════════════════════════════════════════════════════════════════
# Panel a — Population PCA
# ══════════════════════════════════════════════════════════════════════════════
ax_pca = fig.add_subplot(gs[0, 0])

# Reference lines at origin
ax_pca.axhline(0, color="#d0d0d0", lw=0.5, ls="--", zorder=1)
ax_pca.axvline(0, color="#d0d0d0", lw=0.5, ls="--", zorder=1)

# Scatter + label each population
label_offsets = {      # (dx_pt, dy_pt)  fine-tuned to avoid overlap
    "BVI":     (-32,  5),
    "TT":      (  5, -10),
    "Bahamas": (  5,   5),
    "Barbados":(-40,   5),
    "Bermuda": (  5,   5),
    "StLucia": (  5, -10),
}

for pop in pops:
    x, y = pca.loc[pop, "PC1"], pca.loc[pop, "PC2"]
    ax_pca.scatter(
        x, y, s=65, color=COLORS[pop],
        edgecolors="k", linewidths=0.55, zorder=4,
    )
    dx, dy = label_offsets.get(pop, (6, 4))
    ann_kw = dict(
        xy=(x, y), xytext=(dx, dy),
        textcoords="offset points",
        fontsize=6, color="#1a1a1a",
        ha="left" if dx >= 0 else "right",
    )
    if abs(dx) > 15:   # add a subtle leader line for offset labels
        ann_kw["arrowprops"] = dict(
            arrowstyle="-", color="#aaaaaa", lw=0.45, shrinkA=0, shrinkB=2
        )
    ax_pca.annotate(DISPLAY[pop].replace("\n", " "), **ann_kw)

ax_pca.set_xlabel(f"PC1 ({var_exp[0]:.1f}% variance explained)", fontsize=7)
ax_pca.set_ylabel(f"PC2 ({var_exp[1]:.1f}% variance explained)", fontsize=7)
ax_pca.set_title("Population structure", fontsize=7.5, pad=3, fontweight="normal")
ax_pca.spines[["top", "right"]].set_visible(False)
_panel_label(ax_pca, "a")


# ══════════════════════════════════════════════════════════════════════════════
# Panel b — Pairwise FST heat map (UPGMA-ordered)
# ══════════════════════════════════════════════════════════════════════════════
ax_hm = fig.add_subplot(gs[0, 1])

fst_ord = fst.loc[ord_pops, ord_pops].values
n       = len(ord_pops)
xlabels = [DISPLAY[p].replace("\n", " ") for p in ord_pops]

im = ax_hm.imshow(
    fst_ord, cmap="RdYlBu_r", aspect="auto",
    vmin=0.0, vmax=0.030, interpolation="nearest",
)

# Thin white cell borders
for k in range(n + 1):
    ax_hm.axhline(k - 0.5, color="white", lw=0.7, zorder=3)
    ax_hm.axvline(k - 0.5, color="white", lw=0.7, zorder=3)

# Annotate each cell
for i in range(n):
    for j in range(n):
        val = fst_ord[i, j]
        if i == j:
            cell_str, fc = "—", "#888888"
        else:
            cell_str = f"{val:.4f}"
            fc = "white" if val > 0.019 else "#1a1a1a"
        ax_hm.text(
            j, i, cell_str,
            ha="center", va="center",
            fontsize=4.5, color=fc,
        )

ax_hm.set_xticks(range(n))
ax_hm.set_yticks(range(n))
ax_hm.set_xticklabels(xlabels, rotation=35, ha="right", fontsize=5.8)
ax_hm.set_yticklabels(xlabels, fontsize=5.8)
ax_hm.tick_params(axis="both", length=0)

cbar = plt.colorbar(im, ax=ax_hm, fraction=0.044, pad=0.03)
cbar.set_label(r"Pairwise $F_{ST}$ (WC84)", fontsize=6, labelpad=3)
cbar.ax.tick_params(labelsize=5.5)
cbar.outline.set_linewidth(0.4)

ax_hm.set_title(r"Pairwise genetic differentiation ($F_{ST}$)", fontsize=7.5, pad=3,
                fontweight="normal")
_panel_label(ax_hm, "b")


# ══════════════════════════════════════════════════════════════════════════════
# Panel c — UPGMA dendrogram
# ══════════════════════════════════════════════════════════════════════════════
ax_dg = fig.add_subplot(gs[0, 2])

dn = dendrogram(
    Z,
    labels=[DISPLAY[p].replace("\n", " ") for p in pops],
    ax=ax_dg,
    orientation="left",
    leaf_font_size=6.0,
    above_threshold_color="#555555",
    color_threshold=0,                   # all branches same colour
    link_color_func=lambda _: "#444444",
)

# Colour leaf tick labels to match scatter
for tick in ax_dg.get_ymajorticklabels():
    lbl = tick.get_text()
    for pop, col in COLORS.items():
        if DISPLAY[pop].replace("\n", " ") == lbl:
            tick.set_color(col)
            tick.set_fontweight("bold")
            break

ax_dg.set_xlabel(r"$F_{ST}$ distance (UPGMA)", fontsize=7)
ax_dg.set_title("Population clustering", fontsize=7.5, pad=3, fontweight="normal")
ax_dg.xaxis.set_major_formatter(plt.FormatStrFormatter("%.4f"))
ax_dg.spines[["top", "right", "left"]].set_visible(False)
ax_dg.tick_params(axis="y", length=0)
_panel_label(ax_dg, "c")


# ── Save ──────────────────────────────────────────────────────────────────────
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

fig.savefig(OUT_PDF, bbox_inches="tight")
fig.savefig(OUT_PNG, bbox_inches="tight", dpi=300)
plt.close(fig)

print(f"Saved PDF → {OUT_PDF}")
print(f"Saved PNG → {OUT_PNG}")
