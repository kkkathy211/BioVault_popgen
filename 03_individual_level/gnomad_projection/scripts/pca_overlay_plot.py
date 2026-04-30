#!/usr/bin/env python
"""
Overlay study samples on gnomAD HGDP+1kGP reference PCA, colored by ancestry.

Reference PC coords come from `hgdp_tgp_sample_meta.tsv`:
  gnomad_population_inference.pca_scores  -> 16 PCs (same space as the
                                              loadings file used by pca_project.py)
  gnomad_population_inference.pop         -> ancestry label (nfe/afr/eas/...)

Study PC coords come from study_pca_projection.tsv produced by pca_project.py.

Usage:
    python pca_overlay_plot.py <meta_tsv> <study_tsv> <out_png> [pcx pcy]
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# gnomAD population-inference labels (v3.1 set used for the HGDP+1kGP PCA)
POP_COLORS = {
    "afr": "#E69F00",   # African
    "amr": "#56B4E9",   # Admixed American
    "asj": "#CC79A7",   # Ashkenazi Jewish
    "eas": "#009E73",   # East Asian
    "fin": "#0072B2",   # Finnish
    "nfe": "#D55E00",   # Non-Finnish European
    "sas": "#F0E442",   # South Asian
    "mid": "#999999",   # Middle Eastern
    "ami": "#7F00FF",   # Amish
    "oth": "#BBBBBB",   # Other / unassigned
}
POP_LABELS = {
    "afr": "African",
    "amr": "Admixed American",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "fin": "Finnish European",
    "nfe": "Non-Finnish European",
    "sas": "South Asian",
    "mid": "Middle Eastern",
    "ami": "Amish",
    "oth": "Other",
}


def parse_ref(meta_tsv: Path):
    rows = []
    df = pd.read_csv(meta_tsv, sep="\t", dtype=str, keep_default_na=False)
    for s, raw in zip(df["s"], df["gnomad_population_inference"]):
        if not raw or raw == "NA":
            continue
        try:
            j = json.loads(raw)
        except json.JSONDecodeError:
            continue
        scores = j.get("pca_scores")
        pop = j.get("pop") or "oth"
        if not scores:
            continue
        rows.append((s, pop, scores))
    sids = [r[0] for r in rows]
    pops = [r[1] for r in rows]
    pcs = np.array([r[2] for r in rows], dtype=float)
    return sids, pops, pcs


def parse_study(study_tsv: Path):
    df = pd.read_csv(study_tsv, sep="\t")
    # `scores` column is a Hail-exported array literal: "[0.01,-0.02,...]"
    scores = np.array([
        json.loads(s) if isinstance(s, str) else list(s)
        for s in df["scores"]
    ], dtype=float)
    return df["s"].tolist(), scores


def main():
    if len(sys.argv) < 4:
        print("Usage: pca_overlay_plot.py <meta_tsv> <study_tsv> <out_png> [pcx pcy]")
        sys.exit(1)

    meta_tsv = Path(sys.argv[1])
    study_tsv = Path(sys.argv[2])
    out_png = Path(sys.argv[3])
    pcx = int(sys.argv[4]) - 1 if len(sys.argv) > 4 else 0
    pcy = int(sys.argv[5]) - 1 if len(sys.argv) > 5 else 1

    print(f"Loading reference PCs from {meta_tsv}")
    ref_sids, ref_pops, ref_pcs = parse_ref(meta_tsv)
    print(f"  {len(ref_sids)} reference samples, {ref_pcs.shape[1]} PCs")

    print(f"Loading study PCs from {study_tsv}")
    study_sids, study_pcs = parse_study(study_tsv)
    print(f"  {len(study_sids)} study samples, {study_pcs.shape[1]} PCs")

    if ref_pcs.shape[1] != study_pcs.shape[1]:
        print(f"WARNING: ref has {ref_pcs.shape[1]} PCs, study has {study_pcs.shape[1]}. "
              f"Using min.")

    fig, ax = plt.subplots(figsize=(10, 8))

    # reference, one color per pop
    pops_present = sorted(set(ref_pops), key=lambda p: list(POP_COLORS).index(p)
                          if p in POP_COLORS else 99)
    for pop in pops_present:
        mask = np.array([p == pop for p in ref_pops])
        ax.scatter(ref_pcs[mask, pcx], ref_pcs[mask, pcy],
                   s=10, alpha=0.55, linewidth=0,
                   color=POP_COLORS.get(pop, "#BBBBBB"),
                   label=f"{POP_LABELS.get(pop, pop)} (n={mask.sum()})")

    # study samples on top
    ax.scatter(study_pcs[:, pcx], study_pcs[:, pcy],
               s=220, marker="*", c="black",
               edgecolors="gold", linewidth=1.0,
               label=f"Study (n={len(study_sids)})", zorder=5)
    for sid, x, y in zip(study_sids, study_pcs[:, pcx], study_pcs[:, pcy]):
        ax.annotate(str(sid), (x, y), fontsize=7, alpha=0.85,
                    xytext=(5, 5), textcoords="offset points", zorder=6)

    ax.set_xlabel(f"PC{pcx + 1}")
    ax.set_ylabel(f"PC{pcy + 1}")
    ax.set_title("Study samples projected onto gnomAD HGDP+1kGP PC space")
    ax.grid(True, alpha=0.2)
    ax.axhline(0, color="grey", linewidth=0.4)
    ax.axvline(0, color="grey", linewidth=0.4)
    ax.legend(loc="best", fontsize=8, framealpha=0.9)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
