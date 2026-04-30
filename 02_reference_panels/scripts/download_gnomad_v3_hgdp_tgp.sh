#!/usr/bin/env bash
# download_gnomad_v3_hgdp_tgp.sh
# -----------------------------
# Download gnomAD v3.1.2 HGDP+1KGP (TGP) joint-call subset. Contains
# ~4,000 reference individuals from HGDP and the 1000 Genomes Project,
# fully phased with metadata for population/superpopulation labels.
#
# Used by:
#   - 03_individual_level/gnomad_projection  (PCA projection onto HGDP+TGP space)
#
# Files: chr1..22 (~80 GB total)
# Requires: curl

set -euo pipefail

OUTDIR="${1:-./gnomad_v3_hgdp_tgp}"
mkdir -p "$OUTDIR"

BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1.2/vcf/genomes"

for chr in {1..22}; do
  for ext in vcf.bgz vcf.bgz.tbi; do
    f="gnomad.genomes.v3.1.2.hgdp_tgp.chr${chr}.${ext}"
    if [[ -f "$OUTDIR/$f" ]]; then
      echo "Skipping (exists): $f"
      continue
    fi
    echo "Downloading: $f"
    curl -L -o "$OUTDIR/$f" "$BASE/$f"
  done
done

echo "Done. Files in: $OUTDIR"
