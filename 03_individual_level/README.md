# 03 — Individual-Level Analyses

Analyses whose primary input is **per-individual genotype data** (one file per person, GSA microarray format from `01_mock_data_generation/output/`).

## The four analyses

| Sub-analysis | What it does | Reference panel needed |
|---|---|---|
| [`pca_qc/`](pca_qc/) | Quality control + unsupervised PCA on the study samples alone. Sanity check that genotypes are well-formed and samples cluster sensibly. | None |
| [`admixture/`](admixture/) (DELETED!!!)| Estimate ancestry proportions per individual against the 1KGP reference (AFR/EUR/EAS/SAS/AMR), plus local-ancestry inference along the genome. | 1KGP high-coverage |
| [`gnomad_projection/`](gnomad_projection/) | Project study samples into the gnomAD HGDP+TGP PCA space; visualize where they sit relative to global populations. | gnomAD v3.1.2 HGDP+TGP |
| [`sex_biased_admixture/`](sex_biased_admixture/) | Test whether X-chromosome AFR ancestry exceeds autosomal AFR ancestry — the genetic signature of sex-biased colonial admixture. | None |

## Common input

All four read individual genotypes from:

```
../01_mock_data_generation/output/{sample_id}/{sample_id}_X_X_GSAv3-DTC_GRCh38-*.txt
```

If you swap in real GSA data, place it in the same folder structure with the same naming convention and the pipelines will pick it up.

## Independence

These analyses are **independent of each other** — running one does not produce inputs for another. They can be run in any order, in parallel, or skipped individually.
