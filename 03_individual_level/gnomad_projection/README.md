# gnomad_projection — PCA Projection onto gnomAD HGDP+TGP Space

Project study samples into the principal-component space defined by the gnomAD v3.1.2 HGDP+TGP joint-call subset (~4,000 reference individuals from HGDP + 1000 Genomes). This places the study cohort on the same axes as a global reference, so you can read off where each sample lands relative to AFR/AMR/EAS/EUR/SAS clusters.

Also runs ADMIXTURE on the merged cohort across `K = 3..10`.

## Inputs

- Per-individual GSA TXT files at `../../01_mock_data_generation/output/{id}/...txt`.
- gnomAD HGDP+TGP VCFs — download via `../../02_reference_panels/scripts/download_gnomad_v3_hgdp_tgp.sh`.

## Reference / panel files (NOT Uploaded HERE!!!)

```
reference/
├── panel_hgdp_tgp.tsv          # study panel definition
├── hgdp_tgp_sample_meta.tsv    # per-sample metadata (population, superpopulation)
├── hgdp_tgp_samples.txt        # sample list
└── pca_loadings/
    └── gnomad.v3.1.pca_loadings.ht/   # Hail Table of pre-computed PCA loadings (skip in git)
```

## Scripts

```
scripts/
├── convert_ddna_to_plink.py        # GSA TXT → PLINK
├── extract_loadings_variants.py    # Extract SNPs that match the gnomAD PCA loadings
├── pca_project.py                  # Project study samples into gnomAD PC space (Hail)
├── pca_overlay_plot.py             # Continental-scale overlay: study + reference
├── pca_overlay_subpop.py           # Sub-population overlay (CSA, EAS, etc.)
├── visualize_results.py            # ADMIXTURE bar plots
├── run_pipeline_gnomad.sh          # Full pipeline: extract → merge → ADMIXTURE → project
└── run_pca_project.sh              # PCA projection only
```

## Run

```bash
cd scripts
bash run_pipeline_gnomad.sh
# or, projection only (after merge already exists)
bash run_pca_project.sh
```

## Outputs

```
results/
├── admixture/
│   ├── study_admixture.{3..10}.Q   # per-sample ancestry proportions
│   ├── study_admixture.{3..10}.P   # per-allele frequencies
│   ├── admixture_K{3..10}_study.tsv
│   ├── cv_errors.txt               # cross-validation error per K
│   └── admixture_K*.log
├── pca/
│   ├── merged_pca.eigenvec         # PCs on the merged cohort
│   ├── merged_pca.eigenval
│   └── study_pca_coordinates.tsv   # study samples in merged PC space
├── pca_projection/
│   ├── study_pca_projection.tsv    # study samples projected onto gnomAD loadings
│   └── hail.log
├── run_summary/                    # bundled "everything you need" snapshot
└── qc_report.txt

plots/
├── admixture_K{3..10}.png          # ADMIXTURE bar charts per K
├── admixture_cv_error.png          # CV error elbow plot
├── pca_continental.png             # superpopulation-colored PCA
├── pca_subpopulation.png           # finer subpopulation labels
├── pca_overlay_pc1_pc2.png         # study cohort overlaid on reference
├── pca_projection.png              # projection-only view
├── pca_subpop_CSA.png
├── pca_subpop_EAS.png
└── pca_subpop_all.png

logs/                               # multiple timestamped run logs
```

## Notes

- This pipeline depends on Hail for the projection step. Make sure you have a working Hail / Spark setup.
- The `pca_loadings/*.ht/` directory is Hail's table format (a directory full of Parquet files); it's gitignored. Re-extract from gnomAD if needed.
