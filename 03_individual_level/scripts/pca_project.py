#!/usr/bin/env python
"""
Project study samples onto the gnomAD HGDP+1kGP reference PC space using
gnomAD's pre-computed PCA loadings.

Why this is fast:
  - The loadings file is a Hail Table (~hundreds of MB) of variant -> PC
    weights, computed once on the full HGDP+1kGP reference panel.
  - Projection is a sparse matrix multiply: PC_study = (G_study - 2*AF) @ L
  - We never download the per-sample reference VCFs.
  - Total runtime: ~minutes after Hail/Spark warm-up.

Outputs:
  <output_dir>/study_pca_projection.tsv   one row per study sample, PC1..PCn
  <output_dir>/pca_projection.png         scatter of study samples in PC space
                                          overlaid on reference centroids

Usage:
    python pca_project.py <study_plink_prefix> <output_dir>
"""

import sys
from pathlib import Path

import hail as hl


# --- gnomAD PCA loadings location -------------------------------------------
# These loadings were computed on the HGDP+1kGP reference panel as part of
# the gnomAD v3.1 ancestry inference workflow. Verify with:
#   gsutil ls gs://gcp-public-data--gnomad/release/3.1/pca/
#
# If this exact path 404s, also try:
#   gs://gcp-public-data--gnomad/release/3.1.2/pca/gnomad.v3.1.2.hgdp_1kg_subset_pop_pca_loadings.ht
LOADINGS_PATHS = [
    "/Users/spectremac/Desktop/kathyproject/pipeline_gnomad/reference/pca_loadings/gnomad.v3.1.pca_loadings.ht",
    "gs://gcp-public-data--gnomad/release/3.1.2/pca/gnomad.v3.1.2.hgdp_1kg_subset_pop_pca_loadings.ht",
]


def main():
    if len(sys.argv) < 3:
        print("Usage: pca_project.py <study_plink_prefix> <output_dir>")
        sys.exit(1)

    study_prefix = sys.argv[1]
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Hail with local Spark (no cluster needed for ~10 samples) ----------
    hl.init(
        default_reference="GRCh38",
        log=str(output_dir / "hail.log"),
        quiet=True,
        spark_conf={
            "spark.driver.memory": "8g",
            "spark.executor.memory": "8g",
        },
    )

    # --- Study data ---------------------------------------------------------
    print(f"Loading study PLINK from {study_prefix}.{{bed,bim,fam}} ...")
    study_mt = hl.import_plink(
        bed=f"{study_prefix}.bed",
        bim=f"{study_prefix}.bim",
        fam=f"{study_prefix}.fam",
        reference_genome="GRCh38",
        contig_recoding={str(c): f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]},
    )
    n_study = study_mt.count_cols()
    n_var = study_mt.count_rows()
    print(f"  Study: {n_study} samples × {n_var:,} variants")

    # --- Reference loadings -------------------------------------------------
    loadings_ht = None
    for path in LOADINGS_PATHS:
        try:
            print(f"Trying loadings: {path}")
            loadings_ht = hl.read_table(path)
            print(f"  OK -> {path}")
            break
        except Exception as e:
            print(f"  failed: {e}")

    if loadings_ht is None:
        print("ERROR: could not load any of the candidate PCA loadings paths.")
        print("List the bucket and update LOADINGS_PATHS at the top of this file:")
        print("  gsutil ls gs://gcp-public-data--gnomad/release/3.1/pca/")
        print("  gsutil ls gs://gcp-public-data--gnomad/release/3.1.2/pca/")
        sys.exit(1)

    # The loadings table normally has columns: locus, alleles, loadings, pca_af
    print(f"  Loadings columns: {list(loadings_ht.row)}")
    n_load = loadings_ht.count()
    print(f"  Loadings: {n_load:,} variants")

    # --- Project ------------------------------------------------------------
    print("Projecting study samples onto reference PC space...")
    projection_ht = hl.experimental.pc_project(
        study_mt.GT,
        loadings_ht.loadings,
        loadings_ht.pca_af,
    )

    out_tsv = output_dir / "study_pca_projection.tsv"
    print(f"Writing PC coordinates -> {out_tsv}")
    projection_ht.export(str(out_tsv))

    # --- Quick sanity print -------------------------------------------------
    df = projection_ht.to_pandas()
    print(f"\nProjected {len(df)} samples")
    print(df.head().to_string(index=False))

    # --- Optional: simple scatter plot (PC1 vs PC2) -------------------------
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        scores = np.array([list(s) for s in df["scores"]])
        if scores.shape[1] >= 2:
            fig, ax = plt.subplots(figsize=(8, 7))
            ax.scatter(scores[:, 0], scores[:, 1],
                       c="black", marker="*", s=200,
                       edgecolors="gold", linewidth=0.8)
            for i, sid in enumerate(df["s"]):
                ax.annotate(str(sid), (scores[i, 0], scores[i, 1]),
                            fontsize=8, alpha=0.7, xytext=(4, 4),
                            textcoords="offset points")
            ax.set_xlabel("PC1 (gnomAD HGDP+1kGP space)")
            ax.set_ylabel("PC2 (gnomAD HGDP+1kGP space)")
            ax.set_title("Study samples projected onto gnomAD reference PCs")
            ax.grid(True, alpha=0.2)
            ax.axhline(0, color="grey", linewidth=0.5)
            ax.axvline(0, color="grey", linewidth=0.5)
            plot_path = output_dir / "pca_projection.png"
            plt.tight_layout()
            plt.savefig(plot_path, dpi=200, bbox_inches="tight")
            plt.close()
            print(f"\nPlot -> {plot_path}")
            print("Note: this plot shows your samples in absolute PC coords. "
                  "To overlay reference centroids, also load the published "
                  "gnomAD reference scores (separate Hail Table).")
    except Exception as e:
        print(f"(skipped plot: {e})")


if __name__ == "__main__":
    main()
