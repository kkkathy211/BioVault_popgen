"""
04_merge_carib_gnomad.py
========================
Build the master allele-frequency matrix used by all downstream analyses.

Input
-----
  ../../fst_pipeline/data/merged/merged_allele_freq_annotated.tsv
      Caribbean per-island allele frequencies.
  ../data/gnomad_v4_af_per_locus.tsv
      gnomAD v4.1 joint global / AFR / NFE / SAS frequencies, locus-keyed.
  ../../fst_pipeline/data/merged/merged_allele_number.tsv  (optional QC)

Output
------
  ../data/master_af_table.tsv
      locus_key  rsid  BVI  TT  Bahamas  Barbados  Bermuda  StLucia
                       gnomAD_global  gnomAD_AFR  gnomAD_NFE  gnomAD_SAS
      One row per SNP that has gnomAD coverage AND non-missing AF in every
      Caribbean cohort.  This is what we feed downstream.
  ../data/master_af_table_summary.txt
      Quick stats on row counts, dropouts, AF distributions per column.
"""

from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path(__file__).resolve().parents[1]
CARIB = BASE.parent / "fst_pipeline" / "data" / "merged" / "merged_allele_freq_annotated.tsv"
GNOMAD = BASE / "data" / "gnomad_v4_af_per_locus.tsv"
OUT    = BASE / "data" / "master_af_table.tsv"
SUMM   = BASE / "data" / "master_af_table_summary.txt"

ISLANDS = ["BVI", "TT", "Bahamas", "Barbados", "Bermuda", "StLucia"]
GNOMAD_COLS = ["gnomAD_global", "gnomAD_AFR", "gnomAD_NFE", "gnomAD_SAS"]


def main() -> None:
    carib = pd.read_csv(CARIB, sep="\t")
    g     = pd.read_csv(GNOMAD, sep="\t")

    n_carib_in = len(carib)
    df = carib.merge(g, on="locus_key", how="left")

    # require gnomAD global + the three reference pops + complete Caribbean data
    mask_carib = df[ISLANDS].notna().all(axis=1)
    mask_g     = df[GNOMAD_COLS].notna().all(axis=1)
    keep = mask_carib & mask_g
    out = df.loc[keep, ["locus_key", "rsid"] + ISLANDS + GNOMAD_COLS].copy()

    out.to_csv(OUT, sep="\t", index=False, float_format="%.6f")

    # write a small summary
    lines = [
        f"Caribbean SNPs in            : {n_carib_in:,}",
        f"After requiring all islands  : {mask_carib.sum():,}",
        f"After requiring gnomAD ref   : {keep.sum():,}",
        f"Final master_af_table rows   : {len(out):,}",
        "",
        "Allele-frequency distributions (mean ± sd):",
    ]
    for c in ISLANDS + GNOMAD_COLS:
        lines.append(f"  {c:<14s} {out[c].mean():.4f} ± {out[c].std():.4f}")
    SUMM.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n[04] wrote {OUT}")
    print(f"[04] wrote {SUMM}")


if __name__ == "__main__":
    main()
