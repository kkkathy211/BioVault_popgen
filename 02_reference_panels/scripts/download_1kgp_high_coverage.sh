#!/usr/bin/env bash
# download_1kgp_high_coverage.sh
# ------------------------------
# Download the 1000 Genomes Project high-coverage (NYGC, 30x) phased panel
# (3,202 individuals). Hosted on the 1000 Genomes FTP at EBI.
#
# Used by:
#   - 03_individual_level/admixture  (admixture / local-ancestry reference)
#
# Files: chr1..22, X  (~1 TB total, plan accordingly)
# Requires: curl

set -euo pipefail

OUTDIR="${1:-./1kgp_high_coverage}"
mkdir -p "$OUTDIR"

BASE="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"

for chr in {1..22} X; do
  for ext in vcf.gz vcf.gz.tbi; do
    f="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.${ext}"
    if [[ -f "$OUTDIR/$f" ]]; then
      echo "Skipping (exists): $f"
      continue
    fi
    echo "Downloading: $f"
    curl -L -o "$OUTDIR/$f" "$BASE/$f"
  done
done

echo "Done. Files in: $OUTDIR"
