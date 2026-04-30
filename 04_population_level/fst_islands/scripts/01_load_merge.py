"""
Step 1 — Load and merge island allele frequency files.

Merge key: locus_key  (chr-pos-ref-alt — unique per file, identical across all files)
           NOT rsid   (rsid is non-unique: 678k unique vs 1.05M rows due to
                        multi-allelic variants sharing the same rsid)

Outputs
-------
  data/merged/merged_allele_freq.tsv    — locus_key × population  (allele_freq)
  data/merged/merged_allele_number.tsv  — locus_key × population  (allele_number)
  data/merged/snp_overlap_summary.txt   — summary statistics
"""

import logging
from pathlib import Path
import pandas as pd

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[1]
RAW_DIR  = BASE_DIR.parent / "raw_allele_freq_country"
OUT_DIR  = BASE_DIR / "data" / "merged"
LOG_DIR  = BASE_DIR / "logs"

# ── Island label → filename ───────────────────────────────────────────────────
ISLAND_FILES = {
    "BVI"      : "allele_freq_BVI.tsv",
    "TT"       : "allele_freq_TT.tsv",
    "Bahamas"  : "allele_freq_bahamas.tsv",
    "Barbados" : "allele_freq_barbados.tsv",
    "Bermuda"  : "allele_freq_bermuda.tsv",
    "StLucia"  : "allele_freq_stlucia.tsv",
}

# ── Logging ───────────────────────────────────────────────────────────────────
LOG_DIR.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "01_load_merge.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)


def load_island(path: Path, label: str) -> pd.DataFrame:
    """
    Read one island TSV. Index on locus_key (unique per row).
    Rows with allele_number == 0 → freq set to NaN (no genotyping data).
    Returns columns: '<label>_freq', '<label>_n', 'rsid'.
    """
    df = pd.read_csv(path, sep="\t",
                     dtype={"locus_key": str, "rsid": str})
    log.info(f"  {label}: {len(df):,} SNPs  |  "
             f"unique locus_key: {df['locus_key'].nunique():,}  |  "
             f"unique rsid: {df['rsid'].nunique():,}")

    # Treat allele_number == 0 as missing (no sample was genotyped here)
    mask_zero = df["allele_number"] == 0
    log.info(f"    {mask_zero.sum():,} rows with allele_number=0 → NaN")
    df.loc[mask_zero, "allele_freq"]   = float("nan")
    df.loc[mask_zero, "allele_number"] = float("nan")

    df = df.set_index("locus_key")

    # Keep allele_freq, allele_number, and rsid (for downstream annotation)
    out = df[["allele_freq", "allele_number"]].rename(
        columns={"allele_freq":   f"{label}_freq",
                 "allele_number": f"{label}_n"}
    )
    out["rsid"] = df["rsid"]          # carry rsid along for reference
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    log.info("Loading island files …")
    frames = []
    for label, fname in ISLAND_FILES.items():
        frames.append(load_island(RAW_DIR / fname, label))

    # ── Inner join on locus_key ─────────────────────────────────────────────
    # All files share identical locus_key universes so the inner join is a
    # no-op on key coverage. NaN rows are removed in the next step.
    log.info("Merging on locus_key (inner join) …")
    # Use lsuffix/rsuffix to keep only one copy of rsid
    merged = frames[0]
    for f in frames[1:]:
        merged = merged.join(f.drop(columns="rsid"), how="inner")

    freq_cols = [c for c in merged.columns if c.endswith("_freq")]
    n_cols    = [c for c in merged.columns if c.endswith("_n")]

    before = len(merged)
    merged = merged.dropna(subset=freq_cols + n_cols)
    after  = len(merged)
    log.info(f"SNPs after dropping any-population missing: "
             f"{after:,}  (removed {before - after:,})")

    # ── Split and save ────────────────────────────────────────────────────────
    freq_matrix = merged[["rsid"] + freq_cols].copy()
    freq_matrix.columns = ["rsid"] + [c.replace("_freq", "") for c in freq_cols]

    n_matrix = merged[n_cols].copy()
    n_matrix.columns = [c.replace("_n", "") for c in n_cols]

    # Pure frequency matrix without rsid column (used by downstream scripts)
    freq_only = freq_matrix.drop(columns="rsid")

    freq_matrix.to_csv(OUT_DIR / "merged_allele_freq_annotated.tsv", sep="\t")
    freq_only.to_csv(OUT_DIR / "merged_allele_freq.tsv",   sep="\t")
    n_matrix.to_csv (OUT_DIR / "merged_allele_number.tsv", sep="\t")

    log.info(f"Saved freq matrix (annotated) → merged_allele_freq_annotated.tsv  {freq_matrix.shape}")
    log.info(f"Saved freq matrix             → merged_allele_freq.tsv             {freq_only.shape}")
    log.info(f"Saved allele-number matrix    → merged_allele_number.tsv           {n_matrix.shape}")

    # ── Summary ───────────────────────────────────────────────────────────────
    pop_names = list(ISLAND_FILES.keys())
    lines = [
        f"Merge key      : locus_key (chr-pos-ref-alt)",
        f"Populations    : {pop_names}",
        f"Overlapping SNPs (all populations non-missing): {after:,}",
        "",
        "Per-population max allele_number (≈ 2 × max individuals genotyped):",
    ]
    for col in n_matrix.columns:
        lines.append(f"  {col:12s}: {int(n_matrix[col].max())}")

    summary_path = OUT_DIR / "snp_overlap_summary.txt"
    summary_path.write_text("\n".join(lines))
    for ln in lines:
        log.info(ln)

    log.info("Step 1 complete.")


if __name__ == "__main__":
    main()
