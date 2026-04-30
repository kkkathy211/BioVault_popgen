#!/usr/bin/env bash
# run_pipeline.sh — orchestrate the popgen_indepth pipeline end-to-end.
#
# Steps (each is independently re-runnable):
#   02  aggregate per-super-pop allele frequencies from the local
#       gnomAD HGDP+1KGP genotype VCFs (no network)
#   04  merge Caribbean + gnomAD into the master_af_table.tsv
#   05  per-island differential-SNP analysis (vs AFR, vs global) + heatmaps
#   06  AIMs panels (AFR/NFE, AFR/SAS) clustermap + AIMs PCA
#
# Notes
# -----
#   * Old 01_extract_loci_bed.py / 02_fetch_gnomad_v4.sh / 03_build_gnomad_af_table.py
#     are kept on disk for documentation but are no longer called: the local
#     HGDP+TGP reference makes the network slice unnecessary.
#   * The intersection between gnomAD's PCA-loadings panel (~5,808 SNPs) and
#     the Caribbean SNPs is ~5,800 — that is the working set for steps 05/06.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PY="${PY:-python3}"

echo "── 02  compute gnomAD AF locally (HGDP+TGP) ──────"
"$PY" 02_compute_gnomad_af_local.py

echo "── 04  merge Caribbean + gnomAD ──────────────────"
"$PY" 04_merge_carib_gnomad.py

echo "── 05  differential SNPs per island ──────────────"
"$PY" 05_differential_snps_per_island.py

echo "── 06  AIMs heatmaps + PCA ───────────────────────"
"$PY" 06_AIMs_dendrogram.py

echo "Done.  Plots: ../plots   Tables: ../data"
