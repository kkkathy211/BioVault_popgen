"""
Step 1: Parse and merge individual genotype files into a single SNP x individuals matrix.

Input:  One folder per individual under DATA_DIR, each containing one .txt file.
        File format (tab-separated, lines starting with '#' are skipped):
            rsid  chromosome  position  genotype  gs  baf  lrr

Output: data/merged/genotype_matrix_raw.tsv   — rsid x sample_id, genotype strings (e.g. "AA")
        data/merged/snp_info.tsv              — rsid, chromosome, position
"""

import os
import glob
import logging
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR   = Path(__file__).resolve().parents[1]          # pipeline/
DATA_DIR   = BASE_DIR.parent                              # out_copy/ — sample folders live here
OUT_DIR    = BASE_DIR / "data" / "merged"
LOG_DIR    = BASE_DIR / "logs"

SAMPLE_DIRS = [
    d for d in DATA_DIR.iterdir()
    if d.is_dir() and d.name not in {"pipeline", "raw_allele_freq_country"}
]

COLS = ["rsid", "chromosome", "position", "genotype", "gs", "baf", "lrr"]

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
LOG_DIR.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "01_merge.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def find_genotype_file(sample_dir: Path) -> Path:
    """Return the first .txt file found in a sample directory."""
    hits = list(sample_dir.glob("*.txt"))
    if not hits:
        raise FileNotFoundError(f"No .txt file in {sample_dir}")
    return hits[0]


def read_genotype_file(path: Path) -> pd.DataFrame:
    """Read a genotype file, skipping comment lines."""
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=COLS,
        dtype={"rsid": str, "chromosome": str, "position": int,
               "genotype": str, "gs": float, "baf": float, "lrr": float},
    )
    # Drop rows with missing rsid or genotype
    df = df.dropna(subset=["rsid", "genotype"])
    # Normalise genotype: uppercase, strip whitespace
    df["genotype"] = df["genotype"].str.upper().str.strip()
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sample_frames = {}
    snp_info = None

    for sample_dir in sorted(SAMPLE_DIRS):
        sample_id = sample_dir.name
        try:
            gfile = find_genotype_file(sample_dir)
            log.info(f"Reading {sample_id}: {gfile.name}")
            df = read_genotype_file(gfile)

            # Keep SNP metadata from the first sample (same array → same SNPs)
            if snp_info is None:
                snp_info = df[["rsid", "chromosome", "position"]].copy()

            sample_frames[sample_id] = df.set_index("rsid")["genotype"].rename(sample_id)

        except Exception as exc:
            log.warning(f"Skipping {sample_id}: {exc}")

    if not sample_frames:
        log.error("No samples loaded — aborting.")
        return

    log.info(f"Loaded {len(sample_frames)} samples.")

    # Merge: outer join so every rsid present in any sample is kept
    matrix = pd.concat(sample_frames.values(), axis=1, join="outer")
    matrix.index.name = "rsid"

    # Save outputs
    matrix_path = OUT_DIR / "genotype_matrix_raw.tsv"
    matrix.to_csv(matrix_path, sep="\t")
    log.info(f"Saved raw genotype matrix → {matrix_path}  shape={matrix.shape}")

    snp_info_path = OUT_DIR / "snp_info.tsv"
    snp_info.drop_duplicates("rsid").to_csv(snp_info_path, sep="\t", index=False)
    log.info(f"Saved SNP info → {snp_info_path}")


if __name__ == "__main__":
    main()
