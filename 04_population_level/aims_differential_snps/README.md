# aims_differential_snps — AIMs and Per-Island Differential SNPs

Identify ancestry-informative markers (AIMs) and per-island differential SNPs by comparing each Caribbean island's allele frequencies against gnomAD reference populations (AFR, NFE, SAS, plus the global cohort).

Two main outputs:
1. **AIMs** — SNPs whose AF differs maximally between two reference populations (e.g., AFR vs NFE), useful for ancestry estimation.
2. **Differential SNPs per island** — SNPs that are enriched or depleted in a specific island vs gnomAD AFR, or vs the gnomAD global cohort.

## Inputs

```
data/
├── caribbean_loci.bed              # BED intervals (chr, start, end) covering loci of interest
├── caribbean_loci_filter.txt       # Filter list (which loci to include)
└── caribbean_loci_keymap.tsv       # locus_key ↔ chr/pos/rsid mapping
```

Plus gnomAD VCFs (download via `../../02_reference_panels/scripts/download_gnomad_v3_sites.sh`).

## Scripts (run in order)

```
scripts/
├── 01_extract_loci_bed.py            # Build BED from caribbean loci of interest
├── 02_compute_gnomad_af_local.py     # Compute gnomAD AF locally (from downloaded sites VCFs)
├── 02_fetch_gnomad_v4.sh             # Alternative: fetch gnomAD v4 via API
├── 03_build_gnomad_af_table.py       # Assemble master AF table (gnomAD populations × loci)
├── 04_merge_carib_gnomad.py          # Merge Caribbean island AF + gnomAD AF on shared loci
├── 05_differential_snps_per_island.py # Find enriched/depleted SNPs per island
├── 06_AIMs_dendrogram.py             # AIMs selection + clustering visualization
└── run_pipeline.sh                   # Orchestrator
```

Scripts `02_compute_gnomad_af_local.py` and `02_fetch_gnomad_v4.sh` are alternatives — pick one (local is faster if you already downloaded gnomAD, API is faster if you only need a few loci).

## Run

```bash
cd scripts
bash run_pipeline.sh
```

## Outputs

```
results/
├── master_af_table.tsv            # gnomAD populations + Caribbean islands × loci
├── master_af_table_summary.txt    # quick stats
├── gnomad_af_qc.tsv               # QC metrics on gnomAD AF extraction
├── gnomad_v4_af_per_locus.tsv     # per-locus gnomAD v4 AF (if v4 path used)
├── aims/
│   ├── aims_AFR_NFE.tsv           # AIMs distinguishing AFR from NFE
│   ├── aims_AFR_SAS.tsv           # AIMs distinguishing AFR from SAS
│   ├── aims_combined.tsv          # union AIMs set
│   └── aims_combined_pca_loadings.tsv
└── differential_snps/
    ├── {island}_vs_AFR_enriched.tsv      # 6 islands × 4 files each = 24 files
    ├── {island}_vs_AFR_depleted.tsv
    ├── {island}_vs_global_enriched.tsv
    ├── {island}_vs_global_depleted.tsv
    └── all_outliers_long.tsv             # consolidated long-format table

plots/
├── aims_AFR_NFE_clustermap.{png,pdf}      # AFR-vs-NFE AIMs heatmap
├── aims_AFR_SAS_clustermap.{png,pdf}      # AFR-vs-SAS AIMs heatmap
├── aims_combined_pca.{png,pdf}            # PCA on the combined AIMs set
├── diff_snps_heatmap_vs_AFR.{png,pdf}     # per-island differential SNPs vs AFR
└── diff_snps_heatmap_vs_global.{png,pdf}  # per-island differential SNPs vs global

logs/
├── 02_compute_gnomad_af_local.log
└── 02_fetch_gnomad_v4.log
```

## Methods note

A detailed walkthrough of how the AIMs and differential SNPs feed into Figure 2B (the multi-panel ancestry figure) is in [`../../docs/methods_figure2B.md`](../../docs/methods_figure2B.md).
