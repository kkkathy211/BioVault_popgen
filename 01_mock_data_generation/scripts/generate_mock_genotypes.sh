#!/usr/bin/env bash
# generate_mock_genotypes.sh
# ---------------------------
# Generate synthetic GSA microarray genotype files using OpenMined biosynth.
# Each run produces N files at: out/{id}/{id}_X_X_GSAv3-DTC_GRCh38-{month}-{day}-{year}.txt
#
# Requirements: Docker (linux/amd64 platform).
# Usage:        bash generate_mock_genotypes.sh
# Output:       ./out/{id}/...txt        (one folder per synthetic individual)
#               ./out/vcf/...vcf.gz      (after the second step)

set -euo pipefail

# --- Step 1: generate synthetic GSA TXT genotype files ---
# --count 10        : number of synthetic individuals
# --alt-frequency 0.5 : uniform diploid signal (no real population structure)
# --seed 100        : reproducibility
docker run --platform linux/amd64 --rm \
  -v "$(pwd):/work" \
  -w /work \
  ghcr.io/openmined/biosynth:latest \
  synthetic \
    --output "out/{id}/{id}_X_X_GSAv3-DTC_GRCh38-{month}-{day}-{year}.txt" \
    --count 10 \
    --threads 10 \
    --alt-frequency 0.5 \
    --seed 100

# --- Step 2: convert the generated TXT genotypes to VCF ---
docker run --platform linux/amd64 --rm \
  -v "$(pwd):/work" \
  -w /work \
  ghcr.io/openmined/biosynth:latest \
  genotype-to-vcf \
    --input out \
    --outdir out/vcf \
    --gzip

echo "Done. TXT genotypes are in ./out/{id}/  ;  VCFs are in ./out/vcf/"
