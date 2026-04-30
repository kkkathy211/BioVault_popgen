#!/usr/bin/env bash
# download_gnomad_v3_sites.sh
# ---------------------------
# Download gnomAD v3.1.2 sites-only VCFs (allele frequencies for the full
# gnomAD cohort, no individual-level genotypes). Hosted on Google Cloud Storage.
#
# Used by:
#   - 04_population_level/aims_differential_snps  (population-level allele frequency comparisons)
#
# Files: chr1..22, X, Y  (~600 GB total, plan accordingly)
# Requires: curl, ~600 GB free disk

set -euo pipefail

OUTDIR="${1:-./gnomad_v3_sites}"
mkdir -p "$OUTDIR"

BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes"

for chr in {1..22} X Y; do
  for ext in vcf.bgz vcf.bgz.tbi; do
    f="gnomad.genomes.v3.1.2.sites.chr${chr}.${ext}"
    if [[ -f "$OUTDIR/$f" ]]; then
      echo "Skipping (exists): $f"
      continue
    fi
    echo "Downloading: $f"
    curl -L -o "$OUTDIR/$f" "$BASE/$f"
  done
done

echo "Done. Files in: $OUTDIR"
