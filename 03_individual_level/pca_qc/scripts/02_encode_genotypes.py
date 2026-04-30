"""
Step 2: Encode raw string genotypes into numeric dosage (0 / 1 / 2) and
        export a PLINK-compatible PED + MAP file pair.

Encoding logic
--------------
For each SNP we determine the two alleles from the data across all samples,
then count the number of copies of the *minor* allele:

    homozygous reference  →  0
    heterozygous          →  1
    homozygous alternate  →  2
    missing / invalid     →  NaN  (written as NA in the matrix, 0 0 in PED)

Input:  data/merged/genotype_matrix_raw.tsv
        data/merged/snp_info.tsv

Output: data/merged/genotype_matrix_numeric.tsv   — dosage matrix (SNP x sample)
        data/plink/genotypes.ped
        data/plink/genotypes.map
"""

import logging
import numpy as np
import pandas as pd
from collections import Counter
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR  = Path(__file__).resolve().parents[1]
MERGED    = BASE_DIR / "data" / "merged"
PLINK_DIR = BASE_DIR / "data" / "plink"
LOG_DIR   = BASE_DIR / "logs"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
LOG_DIR.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    handlers=[
        logging.FileHandler(LOG_DIR / "02_encode.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
VALID_ALLELES = set("ACGT")


def parse_alleles(gt: str):
    """Return (a1, a2) tuple or (None, None) if invalid/missing."""
    if not isinstance(gt, str) or len(gt) != 2:
        return None, None
    a1, a2 = gt[0], gt[1]
    if a1 not in VALID_ALLELES or a2 not in VALID_ALLELES:
        return None, None
    return a1, a2


def infer_alleles(row: pd.Series):
    """Infer REF/ALT alleles for a SNP row from observed genotypes."""
    allele_counts: Counter = Counter()
    for gt in row:
        a1, a2 = parse_alleles(gt)
        if a1 is not None:
            allele_counts[a1] += 1
            allele_counts[a2] += 1

    if len(allele_counts) < 1:
        return None, None
    sorted_alleles = allele_counts.most_common()  # most common = major
    major = sorted_alleles[0][0]
    minor = sorted_alleles[1][0] if len(sorted_alleles) > 1 else major
    return major, minor          # REF=major, ALT=minor


def encode_genotype(gt: str, ref: str, alt: str) -> float:
    """Return dosage of alt allele: 0, 1, or 2. NaN if missing."""
    a1, a2 = parse_alleles(gt)
    if a1 is None:
        return np.nan
    alleles = [a1, a2]
    unknown = set(alleles) - {ref, alt}
    if unknown:
        return np.nan
    return float(alleles.count(alt))


def gt_to_ped_alleles(gt: str, ref: str, alt: str) -> str:
    """Return PED allele string 'A1 A2'. Missing → '0 0'."""
    a1, a2 = parse_alleles(gt)
    if a1 is None or (set([a1, a2]) - {ref, alt}):
        return "0 0"
    return f"{a1} {a2}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    PLINK_DIR.mkdir(parents=True, exist_ok=True)

    log.info("Loading raw genotype matrix …")
    matrix = pd.read_csv(MERGED / "genotype_matrix_raw.tsv", sep="\t", index_col="rsid")
    snp_info = pd.read_csv(MERGED / "snp_info.tsv", sep="\t",
                           dtype={"chromosome": str}).set_index("rsid")

    samples = list(matrix.columns)
    log.info(f"  {matrix.shape[0]} SNPs × {len(samples)} samples")

    # ------------------------------------------------------------------
    # Numeric encoding
    # ------------------------------------------------------------------
    log.info("Inferring alleles and encoding …")
    allele_map = {}            # rsid → (ref, alt)
    numeric_rows = []

    for rsid, row in matrix.iterrows():
        ref, alt = infer_alleles(row)
        allele_map[rsid] = (ref, alt)
        if ref is None:
            numeric_rows.append([np.nan] * len(samples))
            continue
        numeric_rows.append([encode_genotype(row[s], ref, alt) for s in samples])

    numeric = pd.DataFrame(numeric_rows, index=matrix.index, columns=samples)
    numeric.index.name = "rsid"

    out_num = MERGED / "genotype_matrix_numeric.tsv"
    numeric.to_csv(out_num, sep="\t")
    log.info(f"Saved numeric matrix → {out_num}")

    # ------------------------------------------------------------------
    # PLINK PED file  (one row per individual)
    # ------------------------------------------------------------------
    # PED columns: FID IID PaternalID MaternalID Sex Phenotype  [allele pairs …]
    # We use sample_id for both FID and IID; unknowns = 0.
    log.info("Writing PLINK PED …")
    ped_path = PLINK_DIR / "genotypes.ped"
    snp_list = list(matrix.index)

    with open(ped_path, "w") as fped:
        for sid in samples:
            ref_alt = [allele_map[r] for r in snp_list]
            gt_vals  = [matrix.at[r, sid] for r in snp_list]
            allele_cols = " ".join(
                gt_to_ped_alleles(gt, ra[0] or "0", ra[1] or "0")
                for gt, ra in zip(gt_vals, ref_alt)
            )
            fped.write(f"{sid} {sid} 0 0 0 -9 {allele_cols}\n")

    log.info(f"Saved PED → {ped_path}")

    # ------------------------------------------------------------------
    # PLINK MAP file  (one row per SNP)
    # ------------------------------------------------------------------
    # MAP columns: chromosome  rsid  genetic_distance  position
    log.info("Writing PLINK MAP …")
    map_path = PLINK_DIR / "genotypes.map"
    with open(map_path, "w") as fmap:
        for rsid in snp_list:
            chrom = str(snp_info.at[rsid, "chromosome"]) if rsid in snp_info.index else "0"
            pos   = snp_info.at[rsid, "position"]         if rsid in snp_info.index else 0
            # Numeric chromosome codes expected by PLINK (X→23, Y→24, MT→26)
            chrom_code = (
                chrom.replace("XY", "25").replace("X", "23")
                      .replace("Y", "24").replace("MT", "26")
            )
            fmap.write(f"{chrom_code}\t{rsid}\t0\t{pos}\n")

    log.info(f"Saved MAP → {map_path}")
    log.info("Step 2 complete.")


if __name__ == "__main__":
    main()
