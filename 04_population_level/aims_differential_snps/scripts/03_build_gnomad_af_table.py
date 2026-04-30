"""
03_build_gnomad_af_table.py
===========================
Stitch the per-chromosome gnomAD slices into a single tidy table keyed by
locus_key (<chrom>-<pos>-<ref>-<alt>) so it can be merged with the Caribbean
allele-frequency matrix.

Input
-----
  ../data/gnomad_v4_raw/chr*.tsv   (CHROM POS REF ALT AF AF_afr AF_nfe AF_sas)
  ../data/caribbean_loci_keymap.tsv

Output
------
  ../data/gnomad_v4_af_per_locus.tsv
      locus_key  gnomAD_global  gnomAD_AFR  gnomAD_NFE  gnomAD_SAS
"""

from pathlib import Path
import pandas as pd
import numpy as np

BASE   = Path(__file__).resolve().parents[1]
RAWDIR = BASE / "data" / "gnomad_v4_raw"
KEYMAP = BASE / "data" / "caribbean_loci_keymap.tsv"
OUT    = BASE / "data" / "gnomad_v4_af_per_locus.tsv"

COLS = ["CHROM", "POS", "REF", "ALT", "AF", "AF_afr", "AF_nfe", "AF_sas"]


def load_one(p: Path) -> pd.DataFrame:
    if p.stat().st_size == 0:
        return pd.DataFrame(columns=COLS)
    df = pd.read_csv(p, sep="\t", header=None, names=COLS, na_values=["."])
    return df


def main() -> None:
    keymap = pd.read_csv(KEYMAP, sep="\t")
    parts = []
    for p in sorted(RAWDIR.glob("chr*.tsv")):
        df = load_one(p)
        if not df.empty:
            parts.append(df)
    g = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=COLS)

    # build locus_key (strip 'chr' prefix to match Caribbean side)
    g["chrom"] = g["CHROM"].str.replace(r"^chr", "", regex=True)
    g["locus_key"] = (
        g["chrom"] + "-" + g["POS"].astype(str) + "-" + g["REF"] + "-" + g["ALT"]
    )

    # join via keymap so we keep the canonical locus_key spelling
    g = g[["locus_key", "AF", "AF_afr", "AF_nfe", "AF_sas"]].rename(columns={
        "AF":     "gnomAD_global",
        "AF_afr": "gnomAD_AFR",
        "AF_nfe": "gnomAD_NFE",
        "AF_sas": "gnomAD_SAS",
    })

    # multi-allelic loci can yield multiple rows per (chrom,pos) — keep the row
    # whose ref/alt matches; drop_duplicates on locus_key is sufficient because
    # the key contains ref+alt.
    g = g.drop_duplicates(subset="locus_key", keep="first")

    out = keymap[["locus_key"]].merge(g, on="locus_key", how="left")

    n_total = len(out)
    n_hit   = out["gnomAD_global"].notna().sum()
    print(f"[03] gnomAD coverage: {n_hit:,} / {n_total:,} loci "
          f"({n_hit/n_total*100:.1f}%) have any AF")

    out.to_csv(OUT, sep="\t", index=False, float_format="%.6f")
    print(f"[03] wrote {OUT}")


if __name__ == "__main__":
    main()
