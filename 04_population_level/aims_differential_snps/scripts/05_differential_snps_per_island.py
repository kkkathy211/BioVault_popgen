"""
05_differential_snps_per_island.py
==================================
Find, for each Caribbean island, the SNPs whose allele frequency is most
strongly enriched OR most strongly depleted vs.

  (a) gnomAD AFR    (continental-African background)
  (b) gnomAD global (gnomAD-wide background)

For each (island, reference, direction) we keep the top-N outliers ranked by
ΔAF = AF_island − AF_reference, restricted to common variants in the reference
(AF_ref >= MIN_REF_AF) so we don't surface noisy near-zero sites.

Outputs
-------
  ../data/differential_snps/<island>_vs_<ref>_<top|bot>.tsv
      Per-island per-direction tables.
  ../data/differential_snps/all_outliers_long.tsv
      All four directions × all islands stacked, with rank, ΔAF, AF columns.
  ../plots/diff_snps_heatmap_vs_AFR.pdf / .png
  ../plots/diff_snps_heatmap_vs_global.pdf / .png
      Heatmap rows = union of the top N enriched + N depleted SNPs across
      every island for that reference.  Columns = islands + the four gnomAD
      reference pops.  Cells = allele frequency, with a divergent ΔAF colorbar
      flanking each row block.

The combined heatmap is the headline visual that extends Figure 4C beyond the
PGx panel onto a population-genetics-driven SNP set.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, leaves_list

BASE       = Path(__file__).resolve().parents[1]
MASTER     = BASE / "data" / "master_af_table.tsv"
DIFF_DIR   = BASE / "data" / "differential_snps"
PLOTS_DIR  = BASE / "plots"
DIFF_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

ISLANDS     = ["BVI", "TT", "Bahamas", "Barbados", "Bermuda", "StLucia"]
GNOMAD_REFS = ["gnomAD_global", "gnomAD_AFR", "gnomAD_NFE", "gnomAD_SAS"]
TOP_N       = 25            # per (island, reference, direction)
MIN_REF_AF  = 0.01          # require AF_ref >= 1% so ΔAF isn't dominated by 0→0.05

DISPLAY = {
    "BVI": "BVI", "TT": "Trinidad &\nTobago", "Bahamas": "Bahamas",
    "Barbados": "Barbados", "Bermuda": "Bermuda", "StLucia": "St. Lucia",
    "gnomAD_global": "gnomAD\nglobal", "gnomAD_AFR": "gnomAD\nAFR",
    "gnomAD_NFE": "gnomAD\nNFE", "gnomAD_SAS": "gnomAD\nSAS",
}


# ── Outlier selection ────────────────────────────────────────────────────────
def collect_outliers(df: pd.DataFrame, reference: str) -> pd.DataFrame:
    """For each island, return top-N enriched + top-N depleted vs reference.

    Symmetric MAF mirroring: a SNP with very high AF in the reference is
    information-rich just like a SNP with very low AF, so we apply the
    common-variant gate to min(AF_ref, 1-AF_ref) instead of AF_ref directly.
    """
    rows = []
    af_ref = df[reference]
    common = (af_ref.clip(0, 1).combine(1 - af_ref, min) >= MIN_REF_AF)
    sub = df.loc[common].copy()
    for isl in ISLANDS:
        d = sub[isl] - af_ref.loc[sub.index]
        s = sub.assign(delta=d, island=isl, reference=reference)
        top = s.nlargest(TOP_N, "delta").assign(direction="enriched",
                                                rank=range(1, TOP_N + 1))
        bot = s.nsmallest(TOP_N, "delta").assign(direction="depleted",
                                                 rank=range(1, TOP_N + 1))
        rows.append(top)
        rows.append(bot)
    return pd.concat(rows, ignore_index=True)


def write_per_island_tables(out: pd.DataFrame, ref_short: str) -> None:
    cols = ["locus_key", "rsid", "delta", "rank"] + ISLANDS + GNOMAD_REFS
    for isl in ISLANDS:
        for direc in ("enriched", "depleted"):
            sub = (out.query("island == @isl and direction == @direc")
                      .sort_values("rank")[cols])
            f = DIFF_DIR / f"{isl}_vs_{ref_short}_{direc}.tsv"
            sub.to_csv(f, sep="\t", index=False, float_format="%.6f")


# ── Heatmap ──────────────────────────────────────────────────────────────────
def diverging_cmap():
    return LinearSegmentedColormap.from_list(
        "carib_div", ["#2166AC", "#F7F7F7", "#B2182B"]
    )


def plot_heatmap(out: pd.DataFrame, master: pd.DataFrame,
                 reference: str, ref_short: str) -> None:
    keys = (out.sort_values(["direction", "island", "rank"])["locus_key"]
              .drop_duplicates()
              .tolist())
    M = master.set_index("locus_key").loc[keys, ISLANDS + GNOMAD_REFS]

    # cluster rows by ΔAF-pattern across islands (so visually similar SNPs cluster)
    delta = M[ISLANDS].sub(M[reference], axis=0).to_numpy()
    if len(delta) > 2:
        Z = linkage(delta, method="average", metric="euclidean")
        order = leaves_list(Z)
        M = M.iloc[order]
        keys = [keys[i] for i in order]

    rsid_lookup = master.set_index("locus_key")["rsid"].to_dict()
    row_labels = [rsid_lookup.get(k) or k for k in keys]

    n_rows = len(M)
    fig_h  = max(6, n_rows * 0.13)
    fig, ax = plt.subplots(figsize=(7.5, fig_h))

    cmap = diverging_cmap()
    norm = TwoSlopeNorm(vmin=0, vcenter=0.25, vmax=1.0)
    im = ax.imshow(M.to_numpy(), aspect="auto", cmap=cmap, norm=norm)

    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels([DISPLAY.get(c, c) for c in M.columns], rotation=0,
                       fontsize=6.5)
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels(row_labels, fontsize=4.5)
    ax.tick_params(axis="x", which="both", length=0)
    ax.tick_params(axis="y", which="both", length=0)

    # mark the reference column with a red border
    ref_idx = M.columns.get_loc(reference)
    for spine in ("top", "bottom"):
        ax.axvline(ref_idx - 0.5, color="black", lw=0.4)
        ax.axvline(ref_idx + 0.5, color="black", lw=0.4)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Allele frequency", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    ax.set_title(
        f"Caribbean per-island top differential SNPs vs {DISPLAY[reference].replace(chr(10),' ')}\n"
        f"(top {TOP_N} enriched + top {TOP_N} depleted per island; AF_ref ≥ {MIN_REF_AF*100:.0f}% common)",
        fontsize=8,
    )
    fig.tight_layout()
    pdf = PLOTS_DIR / f"diff_snps_heatmap_vs_{ref_short}.pdf"
    png = PLOTS_DIR / f"diff_snps_heatmap_vs_{ref_short}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[05] wrote {pdf}")
    print(f"[05] wrote {png}")


def main() -> None:
    master = pd.read_csv(MASTER, sep="\t")
    print(f"[05] master table rows: {len(master):,}")

    long_rows = []
    for ref, short in [("gnomAD_AFR", "AFR"), ("gnomAD_global", "global")]:
        print(f"[05] ── differential SNPs vs {ref} ──")
        out = collect_outliers(master, ref)
        write_per_island_tables(out, short)
        long_rows.append(out)
        plot_heatmap(out, master, ref, short)

    long = pd.concat(long_rows, ignore_index=True)
    long_path = DIFF_DIR / "all_outliers_long.tsv"
    long.to_csv(long_path, sep="\t", index=False, float_format="%.6f")
    print(f"[05] wrote {long_path}")


if __name__ == "__main__":
    main()
