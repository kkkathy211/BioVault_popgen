# BioVault_popgen

Population-genetics analyses of Caribbean ancestry, organized by **input**: 
individual-level genotype analyses on one side, 
population-level allele-frequency analyses on the other. 
The repo also documents how the input data are generated and which reference panels are used.

## Directory layout

```
BioVault_popgen/
├── 01_mock_data_generation/   # how to generate synthetic GSA genotypes (biosynth)
├── 02_reference_panels/       # download scripts + .tbi indices for gnomAD / 1KGP
├── 03_individual_level/       # analyses whose input is per-individual genotypes
│   ├── pca_qc/                # genotype QC + Python/PLINK PCA
│   ├── admixture/             # 1KGP-based admixture & local ancestry
│   ├── gnomad_projection/     # PCA projection onto gnomAD HGDP+TGP space
│   └── sex_biased_admixture/  # X-chromosome vs autosome ancestry by sex
├── 04_population_level/       # analyses whose input is per-population allele frequencies
│   ├── fst_islands/           # FST between Caribbean islands
│   └── aims_differential_snps/# ancestry-informative markers + per-island differential SNPs
└── docs/                      # see the shared link I posted in WhatsApp 
```

## Workflow

```
              01_mock_data_generation/output/      02_reference_panels/
              (per-individual GSA TXT files)       (gnomAD, 1KGP — download)
                          │                                │
                ┌─────────┴─────────────┐                  │
                ▼                       ▼                  ▼
        03_individual_level/                     04_population_level/
        - pca_qc                                 - fst_islands
        - admixture        ◄─────────────────────  (uses per-country
        - gnomad_projection ◄────────────────────   allele freq TSVs)
        - sex_biased_admixture                   - aims_differential_snps
                                                   (uses gnomAD AF)
```

Numbers are dependency order, not strict sequence — `03` and `04` are independent and can be run in either order once the inputs in `01`/`02` are in place.

## Quick start

1. **Generate mock genotypes** (or supply your own real GSA files in the same naming scheme):
   ```bash
   cd 01_mock_data_generation/scripts
   bash generate_mock_genotypes.sh
   # outputs land in ../output/{id}/{id}_X_X_GSAv3-DTC_GRCh38-...txt
   ```

2. **Download reference panels** you need (skip the ones you don't):
   ```bash
   cd 02_reference_panels/scripts
   bash download_gnomad_v3_hgdp_tgp.sh   # for gnomad_projection
   bash download_gnomad_v3_sites.sh      # for aims_differential_snps
   bash download_1kgp_high_coverage.sh   # for admixture
   ```

3. **Run any analysis.** Each analysis directory is self-contained with its own `run_pipeline.sh` (or equivalent) and a `README.md` explaining inputs, outputs, and dependencies.

## What each analysis answers

| Analysis | Question | Input granularity |
|---|---|---|
| `pca_qc` | Do the samples cluster sensibly after QC? | Individual genotypes |
| `admixture` | What fraction of each individual's genome comes from AFR / EUR / NAT ancestral populations? | Individual genotypes + 1KGP reference |
| `gnomad_projection` | Where do study samples land in the gnomAD HGDP+TGP PCA space? | Individual genotypes + gnomAD HGDP+TGP |
| `sex_biased_admixture` | Is X-chromosome AFR ancestry > autosomal AFR ancestry (signature of sex-biased colonial admixture)? | Individual genotypes |
| `fst_islands` | How genetically differentiated are the Caribbean islands from each other? | Per-island allele frequencies |
| `aims_differential_snps` | Which SNPs best distinguish Caribbean islands from gnomAD reference populations? | Per-island AF vs gnomAD AF |

## Data not in this repo

The following are **not committed** (see `.gitignore`) but can be regenerated or re-downloaded:

- Synthetic GSA TXT files in `01_mock_data_generation/output/*/` — re-create via the biosynth script.
- gnomAD / 1KGP VCFs in `02_reference_panels/` — re-download via the scripts in `scripts/`.
- Large per-island allele frequency TSVs in `04_population_level/fst_islands/data/` — provided by collaborators.
- PLINK binary files (`*.bed`, `*.bim`, `*.fam`) and Hail tables (`*.ht/`) — produced by the pipelines.

## Repository convention

Each analysis directory follows the same skeleton:

```
<analysis>/
├── README.md     # what the analysis does, how to run, expected I/O
├── scripts/      # numbered scripts (01_*, 02_*, ...) + run_pipeline.sh
├── data/         # primary inputs to this analysis (gitignored if large)
├── results/      # tabular outputs (TSV, PLINK, eigenvecs, etc.)
├── plots/        # figures (PNG, PDF)
└── logs/         # pipeline run logs
```

Not every analysis has every folder 
e.g. `pca_qc` reads its inputs from `01_mock_data_generation/output/`, so it doesn't have its own `data/`.
