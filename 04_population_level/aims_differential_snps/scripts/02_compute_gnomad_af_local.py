"""
02_compute_gnomad_af_local.py
=============================
Aggregate per-super-population allele frequencies from the gnomAD HGDP+1KGP
genotype VCFs that were already sliced to gnomAD's PCA-loadings ancestry-
informative panel (~5,808 SNPs) by the kathyproject_copy/pipeline_gnomad
pipeline.  Sample-level super_pop labels come from panel_hgdp_tgp.tsv.

Super-pop mapping (HGDP+TGP label → gnomAD-style label used downstream):
    AFR  (HGDP+TGP)                    → AFR
    EUR  (HGDP+TGP, pop != "FIN")      → NFE   (Non-Finnish European)
    CSA  (HGDP+TGP)                    → SAS   (Central+South Asian, gnomAD groups them as SAS/CSA)
    all samples in panel               → global

Replaces the old network-tabix steps (01_extract_loci_bed.py + 02_fetch_gnomad_v4.sh + 03_build_gnomad_af_table.py).

Input
-----
  ../../kathyproject_copy/pipeline_gnomad/reference/panel_hgdp_tgp.tsv
  ../../kathyproject_copy/pipeline_gnomad/reference/hgdp_tgp_study_snps_chr{1..22}.vcf.gz

Output
------
  ../data/gnomad_v4_af_per_locus.tsv
      locus_key   gnomAD_global   gnomAD_AFR   gnomAD_NFE   gnomAD_SAS
      (column names kept identical to the previous schema so that 04/05/06
       require no edits.)
  ../data/gnomad_af_qc.tsv
      Per-locus AC/AN/AF for each super-pop (audit trail).
  ../logs/02_compute_gnomad_af_local.log
"""

from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd

BASE   = Path(__file__).resolve().parents[1]
REFDIR = BASE.parent / "kathyproject_copy" / "pipeline_gnomad" / "reference"
PANEL  = REFDIR / "panel_hgdp_tgp.tsv"
OUT    = BASE / "data" / "gnomad_v4_af_per_locus.tsv"
QC     = BASE / "data" / "gnomad_af_qc.tsv"
LOG    = BASE / "logs" / "02_compute_gnomad_af_local.log"

CHROMS = list(range(1, 23))


def log(msg: str) -> None:
    print(msg, flush=True)
    with LOG.open("a") as f:
        f.write(msg + "\n")


def gt_to_alt(arr: np.ndarray) -> np.ndarray:
    """Map a 2-D string array of GT calls to per-cell alt-allele dosage.
    0/0,0|0 → 0;  het (any phase/order) → 1;  1/1,1|1 → 2;  anything else → NaN.
    """
    out = np.full(arr.shape, np.nan, dtype=np.float32)
    out[(arr == "0/0") | (arr == "0|0")] = 0.0
    out[(arr == "0/1") | (arr == "0|1") | (arr == "1/0") | (arr == "1|0")] = 1.0
    out[(arr == "1/1") | (arr == "1|1")] = 2.0
    return out


def parse_chrom_vcf(vcf: Path):
    """Return (df_meta, gts) where df_meta has CHROM/POS/REF/ALT/locus_key
    and gts is an (n_var, n_sample) string ndarray, with the matching sample list."""
    samples = subprocess.check_output(
        ["bcftools", "query", "-l", str(vcf)], text=True
    ).strip().split("\n")
    raw = subprocess.check_output(
        ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", str(vcf)],
        text=True,
    ).strip()
    if not raw:
        return None, None, samples
    rows = [line.split("\t") for line in raw.split("\n")]
    df = pd.DataFrame(rows, columns=["CHROM", "POS", "REF", "ALT"] + samples)
    chrom_bare = df["CHROM"].str.replace("^chr", "", regex=True)
    df["locus_key"] = chrom_bare + "-" + df["POS"] + "-" + df["REF"] + "-" + df["ALT"]
    gts = df[samples].to_numpy()
    return df[["CHROM", "POS", "REF", "ALT", "locus_key"]].copy(), gts, samples


def main() -> None:
    LOG.write_text("")  # truncate
    log(f"[02] reading panel: {PANEL}")
    panel = pd.read_csv(PANEL, sep="\t")
    log(f"[02] panel size: {len(panel):,}")
    log("[02] super_pop counts:\n" + panel["super_pop"].value_counts().to_string())

    groups = {
        "AFR":    set(panel.loc[panel.super_pop == "AFR", "sample"]),
        "NFE":    set(panel.loc[(panel.super_pop == "EUR") & (panel["pop"] != "FIN"), "sample"]),
        "SAS":    set(panel.loc[panel.super_pop == "CSA", "sample"]),
        "global": set(panel["sample"]),
    }
    log("[02] target group sizes (in panel):")
    for g, s in groups.items():
        log(f"        {g:<6s} {len(s):>5d}")

    parts = []
    qc_parts = []

    for chr_ in CHROMS:
        vcf = REFDIR / f"hgdp_tgp_study_snps_chr{chr_}.vcf.gz"
        if not vcf.exists():
            log(f"[02] chr{chr_}: VCF missing, skipping")
            continue
        meta, gts, samples = parse_chrom_vcf(vcf)
        if meta is None:
            log(f"[02] chr{chr_}: empty VCF")
            continue
        # intersect each group with samples actually in this VCF
        sample_idx = {s: i for i, s in enumerate(samples)}
        alt = gt_to_alt(gts)

        out_cols = {"locus_key": meta["locus_key"].values}
        qc_cols = {"locus_key": meta["locus_key"].values}
        for g, samp_set in groups.items():
            idxs = [sample_idx[s] for s in samp_set if s in sample_idx]
            if not idxs:
                af = np.full(len(meta), np.nan, dtype=np.float32)
                ac = np.full(len(meta), 0, dtype=np.int32)
                an = np.full(len(meta), 0, dtype=np.int32)
            else:
                sub = alt[:, idxs]
                with np.errstate(invalid="ignore"):
                    ac = np.nansum(sub, axis=1)
                    an = (~np.isnan(sub)).sum(axis=1) * 2
                    af = np.where(an > 0, ac / an, np.nan).astype(np.float32)
            out_cols[f"gnomAD_{g}"] = af
            qc_cols[f"AC_{g}"] = ac.astype(np.int32)
            qc_cols[f"AN_{g}"] = an.astype(np.int32)
            qc_cols[f"AF_{g}"] = af
        parts.append(pd.DataFrame(out_cols))
        qc_parts.append(pd.DataFrame(qc_cols))
        log(f"[02] chr{chr_}: {len(meta):,} variants processed "
            f"(samples in VCF={len(samples)})")

    final = pd.concat(parts, ignore_index=True).drop_duplicates(
        subset="locus_key", keep="first"
    )
    qc    = pd.concat(qc_parts, ignore_index=True).drop_duplicates(
        subset="locus_key", keep="first"
    )

    # rename to match the schema 04 expects: gnomAD_global / _AFR / _NFE / _SAS
    final = final.rename(columns={
        "gnomAD_global": "gnomAD_global",
        "gnomAD_AFR":    "gnomAD_AFR",
        "gnomAD_NFE":    "gnomAD_NFE",
        "gnomAD_SAS":    "gnomAD_SAS",
    })[["locus_key", "gnomAD_global", "gnomAD_AFR", "gnomAD_NFE", "gnomAD_SAS"]]

    final.to_csv(OUT, sep="\t", index=False, float_format="%.6f")
    qc.to_csv(QC, sep="\t", index=False, float_format="%.6f")

    log(f"[02] wrote {len(final):,} loci → {OUT}")
    log(f"[02] wrote QC table          → {QC}")
    log("[02] AF distributions:")
    for c in ["gnomAD_global", "gnomAD_AFR", "gnomAD_NFE", "gnomAD_SAS"]:
        col = final[c].dropna()
        log(f"        {c:<14s} n={len(col):>5d}  mean={col.mean():.4f}  sd={col.std():.4f}")


if __name__ == "__main__":
    main()
