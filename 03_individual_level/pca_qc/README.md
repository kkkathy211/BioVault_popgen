# pca_qc — Genotype QC and Unsupervised PCA

A lightweight Python + PLINK pipeline to validate that the study genotypes are well-formed and to produce a baseline PCA of the study samples on their own (no external reference).

## Inputs

10 individual GSA TXT files at `../../01_mock_data_generation/output/{id}/{id}_X_X_GSAv3-DTC_GRCh38-*.txt`.

## Scripts (run in order)

```
scripts/
├── 01_merge_genotypes.py    # Merge 10 per-individual TXTs → SNP × individual matrix
├── 02_encode_genotypes.py   # Encode genotypes as numeric (0/1/2/NA)
├── 03_plink_qc_pca.sh       # Convert to PLINK PED/MAP, run PLINK QC + PCA
├── 03b_python_pca.py        # Alternative: pure-Python PCA (no PLINK dependency)
├── 04_plot_pca.py           # Scatter plots: PC1×PC2 and PC3×PC4
└── run_pipeline.sh          # Orchestrates all of the above
```

`03_plink_qc_pca.sh` and `03b_python_pca.py` are two alternative implementations of the same step — pick one. PLINK is faster and standard; the Python version is dependency-free.

## Run

```bash
cd scripts
bash run_pipeline.sh
```

## Outputs

```
results/ (**Did not upload here**)
├── merged/
│   ├── genotype_matrix_raw.tsv      # SNP × sample, raw alleles
│   ├── genotype_matrix_numeric.tsv  # SNP × sample, 0/1/2/NA
│   └── snp_info.tsv                 # rsid, chrom, pos
├── plink/
│   ├── genotypes.ped
│   └── genotypes.map
└── pca/
    ├── pca.eigenvec                 # per-sample PC scores
    ├── pca.eigenval                 # variance explained
    └── sample_metadata.tsv

plots/
├── pca_pc1_pc2.png
└── pca_pc3_pc4.png
```

## Note on the synthetic data

Because the synthetic samples in `01_mock_data_generation` are generated with `--alt-frequency 0.5` (no population structure), the PCA scatter is expected to look like a random cloud — there's no real signal to recover. Use this pipeline as a smoke test against real data.
