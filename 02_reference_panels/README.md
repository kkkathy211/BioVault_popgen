# 02 — Reference Panels

External population-genetics reference data used by the analyses in `03_individual_level/` and `04_population_level/`. We do **not** commit the actual VCFs (each panel is tens to hundreds of GB); only the small `.tbi` indices are kept, and the downloads are scripted.

## Layout

```
02_reference_panels/
├── scripts/
│   ├── download_gnomad_v3_sites.sh        # gnomAD v3.1.2 sites-only (allele freqs)
│   ├── download_gnomad_v3_hgdp_tgp.sh     # gnomAD v3.1.2 HGDP+TGP joint-call subset
│   └── download_1kgp_high_coverage.sh     # 1000 Genomes 30x phased panel (3,202)
└── indices/
    ├── gnomad_v3_sites/         # 24 .tbi  (chr 1-22, X, Y)
    ├── gnomad_v3_hgdp_tgp/      # tbi indices (incomplete; mostly chr22 sample)
    └── 1kgp_high_coverage/      # 1 .tbi   (chr1 sample)
```

## What each panel provides

| Panel | What it is | Used by |
|---|---|---|
| **gnomAD v3.1.2 sites** | Allele frequencies for ~76k whole genomes from gnomAD, broken down by population (AFR, AMR, EAS, NFE, SAS, etc.). No individual-level genotypes. | `04_population_level/aims_differential_snps` |
| **gnomAD v3.1.2 HGDP+TGP** | Joint-called subset (~4,000 individuals) of the Human Genome Diversity Project + 1000 Genomes Project, with population/superpopulation labels. Has individual genotypes. | `03_individual_level/gnomad_projection` |
| **1KGP high-coverage** | 1000 Genomes Project NYGC 30× phased panel — 3,202 individuals across 26 populations. | `03_individual_level/admixture` |

## Downloading

```bash
cd scripts
bash download_gnomad_v3_sites.sh        # ~600 GB
bash download_gnomad_v3_hgdp_tgp.sh     # ~80 GB
bash download_1kgp_high_coverage.sh     # ~1 TB
```

Each script accepts an optional output directory: `bash download_gnomad_v3_sites.sh /path/to/store`.

By default the scripts download into the current directory; point downstream pipelines at the resulting path.

## Why are the indices kept?

`.tbi` files are tiny (~few hundred KB each) and are useful as a manifest — they tell you which chromosomes have already been fetched, and tools like `bcftools` and `tabix` can use them to do random-access reads against locally-mounted or remote-hosted VCFs.
