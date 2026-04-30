# 01 — Mock Data Generation

Synthetic GSA (Global Screening Array) microarray genotypes used as input for every individual-level analysis in this repo.

## Why mock data?

Real direct-to-consumer (DTC) genotype data carries privacy constraints. The pipelines in `03_individual_level/` are developed and tested against synthetic data produced by [OpenMined biosynth](https://github.com/OpenMined/biosynth), which mimics the GSA v3 file format (chromosome, position, alleles, genotype calls) without any real person's DNA.

The synthetic data here is generated with `--alt-frequency 0.5`, which means **the data has no real population structure** — all SNPs are uniformly heterozygous at the population level. This is fine for pipeline development (verifying the code runs end-to-end and produces well-formed outputs) but means that real biological signal (admixture, ancestry differentiation) will only show up when these pipelines are re-run on actual samples.

## How to (re)generate

```bash
cd scripts
bash generate_mock_genotypes.sh
```

Requires Docker (linux/amd64 platform). The script runs two stages:

1. `biosynth synthetic` — generates `out/{id}/{id}_X_X_GSAv3-DTC_GRCh38-{date}.txt` for `--count 10` individuals.
2. `biosynth genotype-to-vcf` — converts the TXT genotypes to VCF.

Defaults: 10 individuals, seed=100, alt-frequency=0.5. Edit the script to change.

## Output

```
output/
├── 305387/305387_X_X_GSAv3-DTC_GRCh38-07-06-2025.txt
├── 351218/...
├── 371396/...
├── 515084/...
├── 525943/...
├── 807622/...
├── 817803/...
├── 821425/...
├── 838338/...
└── 959836/...
```

Each TXT file is ~30 MB. The 10 numeric IDs are the canonical sample IDs referenced by all downstream individual-level analyses.

The TXT files themselves are gitignored — re-run `generate_mock_genotypes.sh` to repopulate.

## Used by

- `03_individual_level/pca_qc`
- `03_individual_level/admixture`
- `03_individual_level/gnomad_projection`
- `03_individual_level/sex_biased_admixture`
