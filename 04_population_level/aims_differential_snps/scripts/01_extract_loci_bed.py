"""
01_extract_loci_bed.py
======================
Build a BED file of Caribbean SNPs that are informative enough to be worth
querying against gnomAD.  We pre-filter to common variants (max-MAF across the
six island cohorts >= MIN_MAX_MAF) since SNPs that are monomorphic in every
island can't be a top "differential" hit anyway and would just slow step 02
down.

Input
-----
  ../../fst_pipeline/data/merged/merged_allele_freq_annotated.tsv

Output
------
  ../data/caribbean_loci.bed              4-col BED: chr<N>\tstart\tend\tlocus_key
  ../data/caribbean_loci_keymap.tsv       locus_key chrom_chr pos ref alt rsid
  ../data/caribbean_loci_filter.txt       small summary of the filter step
"""

from pathlib import Path
import pandas as pd
import numpy as np

BASE   = Path(__file__).resolve().parents[1]
SRC    = BASE.parent / "fst_pipeline" / "data" / "merged" / "merged_allele_freq_annotated.tsv"
OUTBED = BASE / "data" / "caribbean_loci.bed"
OUTMAP = BASE / "data" / "caribbean_loci_keymap.tsv"
OUTSUM = BASE / "data" / "caribbean_loci_filter.txt"

ISLANDS = ["BVI", "TT", "Bahamas", "Barbados", "Bermuda", "StLucia"]
CHR_ORDER = {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24, "MT": 25}
MIN_MAX_MAF = 0.05   # keep SNP if at least one island has MAF >= 5%


def main() -> None:
    df = pd.read_csv(SRC, sep="\t")
    n_in = len(df)

    # max MAF across islands
    af = df[ISLANDS].clip(0, 1)
    maf = pd.concat([af, 1 - af]).groupby(level=0).min()  # element-wise min(AF, 1-AF)
    df["max_maf"] = maf.max(axis=1)
    df_keep = df[df["max_maf"] >= MIN_MAX_MAF].copy()

    parts = df_keep["locus_key"].str.split("-", n=3, expand=True)
    parts.columns = ["chrom", "pos", "ref", "alt"]
    parts["pos"] = parts["pos"].astype(int)
    df_keep = pd.concat([df_keep[["locus_key", "rsid"]].reset_index(drop=True),
                         parts.reset_index(drop=True)], axis=1)

    df_keep = df_keep[df_keep["chrom"].isin(CHR_ORDER)].copy()
    df_keep["chrom_chr"]  = "chr" + df_keep["chrom"]
    df_keep["chrom_sort"] = df_keep["chrom"].map(CHR_ORDER)
    df_keep = (df_keep.sort_values(["chrom_sort", "pos"])
                       .drop(columns="chrom_sort"))

    bed = df_keep[["chrom_chr", "pos", "pos", "locus_key"]].copy()
    bed.columns = ["chrom", "start", "end", "name"]
    bed["start"] = bed["start"] - 1
    bed.to_csv(OUTBED, sep="\t", header=False, index=False)

    keymap = df_keep[["locus_key", "chrom_chr", "pos", "ref", "alt", "rsid"]]
    keymap.to_csv(OUTMAP, sep="\t", index=False)

    n_out = len(bed)
    summary = (
        f"Input SNPs                       : {n_in:,}\n"
        f"MAF filter threshold (max across islands): {MIN_MAX_MAF}\n"
        f"SNPs passing filter              : {n_out:,}  "
        f"({n_out/n_in*100:.1f}%)\n"
        f"BED file                         : {OUTBED}\n"
        f"Keymap                           : {OUTMAP}\n"
    )
    OUTSUM.write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
