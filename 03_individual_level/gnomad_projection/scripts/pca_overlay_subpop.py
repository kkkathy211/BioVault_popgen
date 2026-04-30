#!/usr/bin/env python
"""
Overlay study samples on gnomAD HGDP+1kGP reference PCA, colored by
fine-grained sub-population (e.g. Han, Japanese, Yoruba, Mbuti, ...).

PC coords come from `gnomad_population_inference.pca_scores` (same space
as study projection). Labels come from `hgdp_tgp_meta`:
  genetic_region  -> AFR/AMR/CSA/EAS/EUR/MID/OCE
  population      -> sub-population name (~80 categories)

Usage:
    python pca_overlay_subpop.py <meta_tsv> <study_tsv> <out_png> \\
           [pcx pcy] [region_filter]

    region_filter: optional one of AFR AMR CSA EAS EUR MID OCE
                   (filters reference to that region only -- useful for
                    zooming into e.g. East Asian sub-pops without 80-color
                    clutter)
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

# Continental region -> matplotlib colormap (palette per region so subpops
# within a region share a hue family).
REGION_CMAPS = {
    "AFR": "Oranges",
    "AMR": "Blues",
    "CSA": "YlGn",
    "EAS": "Greens",
    "EUR": "Reds",
    "MID": "Purples",
    "OCE": "PuRd",
}
REGION_ORDER = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID", "OCE"]


def parse_ref(meta_tsv: Path):
    df = pd.read_csv(meta_tsv, sep="\t", dtype=str, keep_default_na=False)
    rows = []
    for s, pop_inf, hgdp in zip(df["s"],
                                 df["gnomad_population_inference"],
                                 df["hgdp_tgp_meta"]):
        if not pop_inf or pop_inf == "NA" or not hgdp or hgdp == "NA":
            continue
        try:
            pi = json.loads(pop_inf)
            hm = json.loads(hgdp)
        except json.JSONDecodeError:
            continue
        scores = pi.get("pca_scores")
        region = hm.get("genetic_region")
        subpop = hm.get("population")
        if not scores or not region or not subpop:
            continue
        rows.append((s, region, subpop, scores))
    return rows


def parse_study(study_tsv: Path):
    df = pd.read_csv(study_tsv, sep="\t")
    scores = np.array([
        json.loads(s) if isinstance(s, str) else list(s)
        for s in df["scores"]
    ], dtype=float)
    return df["s"].tolist(), scores


def assign_subpop_colors(rows, region_filter=None):
    """Group subpops by region; assign a shade within the region's colormap."""
    by_region = {}
    for _, region, subpop, _ in rows:
        if region_filter and region != region_filter:
            continue
        by_region.setdefault(region, set()).add(subpop)

    colors = {}
    for region in REGION_ORDER:
        if region not in by_region:
            continue
        subpops = sorted(by_region[region])
        cmap = cm.get_cmap(REGION_CMAPS.get(region, "Greys"))
        # Sample within 0.35..0.95 to avoid invisible-pale and pure-black.
        for i, sp in enumerate(subpops):
            shade = 0.35 + (0.6 * (i / max(1, len(subpops) - 1)))
            colors[(region, sp)] = cmap(shade)
    return colors, by_region


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    meta_tsv = Path(sys.argv[1])
    study_tsv = Path(sys.argv[2])
    out_png = Path(sys.argv[3])
    pcx = int(sys.argv[4]) - 1 if len(sys.argv) > 4 else 0
    pcy = int(sys.argv[5]) - 1 if len(sys.argv) > 5 else 1
    region_filter = sys.argv[6].upper() if len(sys.argv) > 6 else None

    print(f"Loading reference from {meta_tsv}")
    rows = parse_ref(meta_tsv)
    print(f"  {len(rows)} HGDP+1kGP reference samples with subpop labels")

    if region_filter:
        rows_plot = [r for r in rows if r[1] == region_filter]
        print(f"  filtered to region={region_filter}: {len(rows_plot)} samples")
    else:
        rows_plot = rows

    print(f"Loading study from {study_tsv}")
    study_sids, study_pcs = parse_study(study_tsv)
    print(f"  {len(study_sids)} study samples")

    colors, by_region = assign_subpop_colors(rows_plot, region_filter)

    fig, ax = plt.subplots(figsize=(12, 9))

    # Plot reference, sorted so that subpops within a region are grouped in legend.
    for region in REGION_ORDER:
        if region not in by_region:
            continue
        for subpop in sorted(by_region[region]):
            mask_xy = [(r[3][pcx], r[3][pcy])
                       for r in rows_plot
                       if r[1] == region and r[2] == subpop]
            xs = [p[0] for p in mask_xy]
            ys = [p[1] for p in mask_xy]
            ax.scatter(xs, ys, s=14, alpha=0.75, linewidth=0,
                       color=colors[(region, subpop)],
                       label=f"[{region}] {subpop} (n={len(xs)})")

    # Study samples on top
    ax.scatter(study_pcs[:, pcx], study_pcs[:, pcy],
               s=240, marker="*", c="black",
               edgecolors="gold", linewidth=1.2,
               label=f"Study (n={len(study_sids)})", zorder=5)
    for sid, x, y in zip(study_sids, study_pcs[:, pcx], study_pcs[:, pcy]):
        ax.annotate(str(sid), (x, y), fontsize=7, alpha=0.85,
                    xytext=(5, 5), textcoords="offset points", zorder=6)

    title = "Study samples projected onto gnomAD HGDP+1kGP PC space"
    if region_filter:
        title += f"  —  zoomed to {region_filter}"
    ax.set_xlabel(f"PC{pcx + 1}")
    ax.set_ylabel(f"PC{pcy + 1}")
    ax.set_title(title)
    ax.grid(True, alpha=0.2)
    ax.axhline(0, color="grey", linewidth=0.4)
    ax.axvline(0, color="grey", linewidth=0.4)

    # Legend off to the right -- can be very tall when not filtered
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=6 if not region_filter else 8,
              ncol=2 if not region_filter else 1,
              frameon=False)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
