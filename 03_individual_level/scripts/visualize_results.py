#!/usr/bin/env python3
"""
Visualize global ancestry results when the reference panel is the
gnomAD v4 HGDP+1kGP subset.

Differences from the 1KG-only version:
  - 7 continental regions (AFR, AMR, CSA, EAS, EUR, MID, OCE) instead of 5
    (CSA = Central/South Asian, MID = Middle Eastern, OCE = Oceanian).
  - ~80 sub-populations with ethnographic labels (YRI, Sardinian, Karitiana, ...).
  - Panel file format is `panel_hgdp_tgp.tsv`: sample\tpop\tsuper_pop\tproject

Outputs:
  1. PCA plots (PC1×PC2, PC1×PC3) coloured by continental region with study
     samples overlaid.
  2. PCA sub-population plot (PC1×PC2) coloured by ~80 sub-pops — the main
     selling point of swapping to HGDP+1kGP.
  3. ADMIXTURE bar plots for each K, with reference samples grouped by region
     and study samples called out at the right edge.
  4. Per-sample ancestry table.

Usage:
    python visualize_results.py <results_dir> <reference_dir>
"""

import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# ---- Continental region colours (gnomAD HGDP+1kGP scheme) ----
REGION_COLORS = {
    "AFR":   "#E41A1C",  # African
    "AMR":   "#FF7F00",  # Admixed American (incl. HGDP Amerindian)
    "CSA":   "#984EA3",  # Central + South Asian (gnomAD merges these)
    "EAS":   "#4DAF4A",  # East Asian
    "EUR":   "#377EB8",  # European
    "MID":   "#A65628",  # Middle Eastern
    "OCE":   "#F781BF",  # Oceanian
    "Study": "#000000",
}

REGION_LABELS = {
    "AFR":   "African",
    "AMR":   "Admixed American",
    "CSA":   "Central/South Asian",
    "EAS":   "East Asian",
    "EUR":   "European",
    "MID":   "Middle Eastern",
    "OCE":   "Oceanian",
    "Study": "Study Samples",
}

# Some 1KG super_pop codes you may still see in older metadata files.
# Map them onto the gnomAD region codes.
LEGACY_REGION_ALIASES = {
    "SAS": "CSA",          # 1KG South Asian -> gnomAD CSA bucket
    "AMR_HGDP": "AMR",
}


def normalise_region(code):
    if not code:
        return None
    code = code.strip()
    return LEGACY_REGION_ALIASES.get(code, code)


# ---- Helpers ----
def read_eigenvec(filepath):
    """plink2 eigenvec -> dict IID -> [PC1, PC2, ...]"""
    data = {}
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if parts[0] in ("#FID", "FID", "#IID"):
                continue
            iid = parts[1]
            try:
                pcs = [float(x) for x in parts[2:]]
            except ValueError:
                continue
            data[iid] = pcs
    return data


def read_panel(reference_dir):
    """Read the HGDP+1kGP panel file produced by the pipeline.
    Returns dict IID -> (pop, super_pop, project)."""
    panel = {}
    panel_file = os.path.join(reference_dir, "panel_hgdp_tgp.tsv")
    if not os.path.exists(panel_file):
        # Fallback: try the raw gnomAD metadata TSV directly.
        for cand in ("hgdp_tgp_sample_meta.tsv", "gnomad_meta_v1.tsv"):
            p = os.path.join(reference_dir, cand)
            if os.path.exists(p):
                panel_file = p
                break

    if not os.path.exists(panel_file):
        print(f"WARNING: no panel file found in {reference_dir}")
        return panel

    with open(panel_file) as f:
        header = f.readline().strip().split("\t")

        def find(*cands):
            for c in cands:
                if c in header:
                    return header.index(c)
            return None

        i_id  = find("sample", "s", "IID", "sample_id")
        i_pop = find("pop", "population", "hgdp_tgp_meta.population")
        i_sup = find("super_pop", "hgdp_tgp_meta.genetic_region",
                     "genetic_region", "superpop")
        i_proj = find("project", "hgdp_tgp_meta.project")

        if i_id is None or i_pop is None or i_sup is None:
            print(f"WARNING: panel file {panel_file} missing expected columns; got {header}")
            return panel

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_id, i_pop, i_sup):
                continue
            iid = parts[i_id]
            pop = parts[i_pop]
            sup = normalise_region(parts[i_sup])
            proj = parts[i_proj] if (i_proj is not None and i_proj < len(parts)) else ""
            panel[iid] = (pop, sup, proj)

    print(f"Loaded panel with {len(panel)} samples from {panel_file}")
    return panel


def read_study_samples(fam_file):
    samples = set()
    if not os.path.exists(fam_file):
        return samples
    with open(fam_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                samples.add(parts[1])
    return samples


def read_eigenvalues(eigenval_file):
    if not os.path.exists(eigenval_file):
        return None
    with open(eigenval_file) as f:
        vals = [float(line.strip()) for line in f if line.strip()]
    total = sum(vals)
    return [100 * v / total for v in vals]


# ---- PCA: continental view ----
def plot_pca_continental(eigenvec, panel, study_ids, var_explained, output_dir):
    region_data = {r: {"pc1": [], "pc2": [], "pc3": []} for r in REGION_COLORS}

    for iid, pcs in eigenvec.items():
        if len(pcs) < 3:
            continue
        if iid in study_ids:
            grp = "Study"
        elif iid in panel:
            _, sup, _ = panel[iid]
            grp = sup
        else:
            continue
        if grp in region_data:
            region_data[grp]["pc1"].append(pcs[0])
            region_data[grp]["pc2"].append(pcs[1])
            region_data[grp]["pc3"].append(pcs[2])

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    for ax_idx, (key, label, idx) in enumerate(
        [("pc2", "PC2", 1), ("pc3", "PC3", 2)]
    ):
        ax = axes[ax_idx]
        for r in ["AFR", "AMR", "CSA", "EAS", "EUR", "MID", "OCE"]:
            d = region_data[r]
            if not d["pc1"]:
                continue
            ax.scatter(
                d["pc1"], d[key],
                c=REGION_COLORS[r],
                label=f'{r} ({REGION_LABELS[r]}, n={len(d["pc1"])})',
                alpha=0.35, s=10, rasterized=True,
            )
        d = region_data["Study"]
        if d["pc1"]:
            ax.scatter(
                d["pc1"], d[key],
                c=REGION_COLORS["Study"],
                label=f'Study (n={len(d["pc1"])})',
                s=130, marker="*",
                edgecolors="gold", linewidth=0.8, zorder=10,
            )
        x_lab = f"PC1 ({var_explained[0]:.1f}%)" if var_explained else "PC1"
        y_lab = f"{label} ({var_explained[idx]:.1f}%)" if var_explained else label
        ax.set_xlabel(x_lab, fontsize=12)
        ax.set_ylabel(y_lab, fontsize=12)
        ax.set_title(f"PCA: PC1 vs {label} (continental)", fontsize=13)
        ax.legend(loc="best", fontsize=8, markerscale=1.5)
        ax.grid(True, alpha=0.2)

    plt.suptitle("PCA — Study Samples vs gnomAD HGDP+1kGP (continental regions)",
                 fontsize=15, y=1.02)
    plt.tight_layout()
    out = os.path.join(output_dir, "pca_continental.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"PCA continental plot: {out}")


# ---- PCA: sub-population view (the gnomAD-specific one) ----
def plot_pca_subpop(eigenvec, panel, study_ids, var_explained, output_dir):
    """Colour by sub-population. Uses one marker per region to keep ~80 pops legible."""
    pop_pts = {}  # pop -> {pc1: [], pc2: [], region: str}
    for iid, pcs in eigenvec.items():
        if len(pcs) < 2 or iid in study_ids or iid not in panel:
            continue
        pop, sup, _ = panel[iid]
        d = pop_pts.setdefault(pop, {"pc1": [], "pc2": [], "region": sup})
        d["pc1"].append(pcs[0])
        d["pc2"].append(pcs[1])

    fig, ax = plt.subplots(figsize=(15, 11))

    # Group sub-pops by region; assign each a colour shade within the region's hue.
    region_to_pops = {}
    for pop, d in pop_pts.items():
        region_to_pops.setdefault(d["region"], []).append(pop)

    markers = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "h", "*"]
    for region, pops in region_to_pops.items():
        base = REGION_COLORS.get(region, "#888888")
        cmap = _shade_palette(base, len(pops))
        for i, pop in enumerate(sorted(pops)):
            d = pop_pts[pop]
            ax.scatter(
                d["pc1"], d["pc2"],
                color=cmap[i],
                marker=markers[i % len(markers)],
                s=18, alpha=0.75, linewidth=0,
                label=f"{pop} ({region}, n={len(d['pc1'])})",
            )

    # Study samples
    sx, sy = [], []
    for iid in study_ids:
        if iid in eigenvec and len(eigenvec[iid]) >= 2:
            sx.append(eigenvec[iid][0])
            sy.append(eigenvec[iid][1])
    if sx:
        ax.scatter(sx, sy, c="#000000", marker="*", s=200,
                   edgecolors="gold", linewidth=0.8,
                   label=f"Study (n={len(sx)})", zorder=20)

    x_lab = f"PC1 ({var_explained[0]:.1f}%)" if var_explained else "PC1"
    y_lab = f"PC2 ({var_explained[1]:.1f}%)" if var_explained else "PC2"
    ax.set_xlabel(x_lab, fontsize=12)
    ax.set_ylabel(y_lab, fontsize=12)
    ax.set_title("PCA by sub-population — gnomAD HGDP+1kGP (~80 populations)",
                 fontsize=14)
    ax.grid(True, alpha=0.2)

    # Legend goes outside; many entries.
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5),
              fontsize=6, ncol=2, frameon=False)

    plt.tight_layout()
    out = os.path.join(output_dir, "pca_subpopulation.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"PCA sub-pop plot: {out}")


def _shade_palette(hex_color, n):
    """Make n shades around a base hex colour (lighter -> darker)."""
    if n <= 1:
        return [hex_color]
    rgb = np.array(matplotlib.colors.to_rgb(hex_color))
    out = []
    for i in range(n):
        # interpolate between 0.45*base+0.55*white and the base colour
        t = 0.4 + 0.6 * (i / max(1, n - 1))
        c = t * rgb + (1 - t) * np.array([1.0, 1.0, 1.0])
        out.append(matplotlib.colors.to_hex(np.clip(c, 0, 1)))
    return out


# ---- ADMIXTURE ----
def plot_admixture(results_dir, reference_dir, output_dir, panel, study_ids):
    admix_dir = os.path.join(results_dir, "admixture")
    if not os.path.isdir(admix_dir):
        print("ADMIXTURE results directory not found")
        return

    q_files = sorted(glob.glob(os.path.join(admix_dir, "*.Q")))
    if not q_files:
        print("No .Q files found")
        return

    fam_file = os.path.join(admix_dir, "study_admixture.fam")
    if not os.path.exists(fam_file):
        print("ADMIXTURE fam file missing")
        return

    with open(fam_file) as f:
        sample_ids = [line.strip().split()[1] for line in f]

    # Region order for grouping the bar plot
    region_order = ["AFR", "MID", "EUR", "CSA", "EAS", "OCE", "AMR"]
    region_idx = {r: i for i, r in enumerate(region_order)}

    def sort_key(i):
        sid = sample_ids[i]
        if sid in study_ids:
            return (len(region_order) + 1, "", sid)  # Study at far right
        if sid in panel:
            pop, sup, _ = panel[sid]
            return (region_idx.get(sup, len(region_order)), pop, sid)
        return (len(region_order), "", sid)

    sort_idx = sorted(range(len(sample_ids)), key=sort_key)
    ids_sorted = [sample_ids[i] for i in sort_idx]

    # CV errors
    cv_errors = {}
    cv_file = os.path.join(admix_dir, "cv_errors.txt")
    if os.path.exists(cv_file):
        with open(cv_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    try:
                        k = int(parts[0].replace("K=", "").replace(":", ""))
                        cv_errors[k] = float(parts[1])
                    except ValueError:
                        pass

    for q_file in q_files:
        try:
            K = int(os.path.basename(q_file).split(".")[-2])
        except (ValueError, IndexError):
            continue

        Q = np.loadtxt(q_file)
        if Q.ndim == 1:
            Q = Q.reshape(-1, 1)
        if Q.shape[0] != len(sample_ids):
            print(f"WARN K={K}: Q rows {Q.shape[0]} != samples {len(sample_ids)}")
            continue

        Q_sorted = Q[sort_idx]

        fig, ax = plt.subplots(figsize=(max(14, len(Q) * 0.012), 5))
        cmap = plt.cm.tab10 if K <= 10 else plt.cm.tab20
        colors = [cmap(i % cmap.N) for i in range(K)]

        bottom = np.zeros(len(Q_sorted))
        for k in range(K):
            ax.bar(range(len(Q_sorted)), Q_sorted[:, k],
                   bottom=bottom, color=colors[k],
                   width=1.0, edgecolor="none",
                   label=f"K{k+1}")
            bottom += Q_sorted[:, k]

        # Region/study separators on the x axis
        prev_region = None
        for x, sid in enumerate(ids_sorted):
            if sid in study_ids:
                region = "Study"
            elif sid in panel:
                _, region, _ = panel[sid]
            else:
                region = "?"
            if region != prev_region:
                ax.axvline(x - 0.5, color="white", linewidth=0.8, alpha=0.7)
                ax.text(x, 1.02, region or "?", fontsize=8, rotation=0,
                        ha="left", va="bottom",
                        color=REGION_COLORS.get(region, "#333"))
                prev_region = region

        ax.set_xlim(-0.5, len(Q_sorted) - 0.5)
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_ylabel("Ancestry proportion", fontsize=12)
        title = f"ADMIXTURE K={K}"
        if K in cv_errors:
            title += f"  (CV error = {cv_errors[K]:.4f})"
        ax.set_title(title, fontsize=13, pad=18)
        ax.legend(loc="upper right", fontsize=7, ncol=min(K, 6))

        plt.tight_layout()
        out = os.path.join(output_dir, f"admixture_K{K}.png")
        plt.savefig(out, dpi=200, bbox_inches="tight")
        plt.close()
        print(f"ADMIXTURE K={K} plot: {out}")

        # Per-study-sample table
        if study_ids:
            tbl_path = os.path.join(output_dir, f"admixture_K{K}_study.tsv")
            with open(tbl_path, "w") as f:
                f.write("SampleID\t" + "\t".join(f"K{i+1}" for i in range(K)) + "\n")
                for sid in sorted(study_ids):
                    if sid in sample_ids:
                        row = Q[sample_ids.index(sid)]
                        f.write(sid + "\t" + "\t".join(f"{v:.4f}" for v in row) + "\n")

    if cv_errors:
        fig, ax = plt.subplots(figsize=(8, 5))
        ks = sorted(cv_errors.keys())
        errs = [cv_errors[k] for k in ks]
        ax.plot(ks, errs, "o-", color="#E41A1C", markersize=8, linewidth=2)
        best_k = ks[int(np.argmin(errs))]
        ax.axvline(best_k, color="gray", linestyle="--", alpha=0.5)
        ax.annotate(f"Best K={best_k}", xy=(best_k, min(errs)),
                    xytext=(best_k + 0.3, min(errs) + 0.001),
                    fontsize=11, color="gray")
        ax.set_xlabel("K", fontsize=12)
        ax.set_ylabel("CV error", fontsize=12)
        ax.set_title("ADMIXTURE Cross-Validation", fontsize=13)
        ax.set_xticks(ks)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        out = os.path.join(output_dir, "admixture_cv_error.png")
        plt.savefig(out, dpi=200, bbox_inches="tight")
        plt.close()
        print(f"CV error plot: {out}")


# ---- Per-study-sample summary ----
def write_study_pca_table(eigenvec, study_ids, output_dir):
    out = os.path.join(output_dir, "study_pca_coordinates.tsv")
    with open(out, "w") as f:
        f.write("SampleID\tPC1\tPC2\tPC3\tPC4\tPC5\n")
        for sid in sorted(study_ids):
            if sid in eigenvec:
                pcs = eigenvec[sid]
                f.write(sid + "\t" + "\t".join(f"{p:.6f}" for p in pcs[:5]) + "\n")
    print(f"Study PCA coords: {out}")


def main():
    if len(sys.argv) < 3:
        print("Usage: visualize_results.py <results_dir> <reference_dir>")
        sys.exit(1)

    results_dir = sys.argv[1]
    reference_dir = sys.argv[2]
    output_dir = os.path.join(results_dir, "plots")
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("ANCESTRY VISUALIZATION — gnomAD HGDP+1kGP reference")
    print("=" * 60)

    eigenvec_file = os.path.join(results_dir, "pca", "merged_pca.eigenvec")
    if not os.path.exists(eigenvec_file):
        print(f"PCA eigenvec missing: {eigenvec_file}")
        sys.exit(0)

    eigenvec = read_eigenvec(eigenvec_file)
    var_explained = read_eigenvalues(
        os.path.join(results_dir, "pca", "merged_pca.eigenval")
    )
    panel = read_panel(reference_dir)

    # Find study samples from the working/ directory
    candidates = [
        os.path.join(results_dir, "..", "working", "study_qc.fam"),
        os.path.join(results_dir, "..", "working", "study_raw.fam"),
    ]
    study_fam = next((p for p in candidates if os.path.exists(p)), None)
    study_ids = read_study_samples(study_fam) if study_fam else set()
    print(f"Study samples: {len(study_ids)}")

    print("\n--- PCA continental ---")
    plot_pca_continental(eigenvec, panel, study_ids, var_explained, output_dir)

    print("\n--- PCA sub-population ---")
    plot_pca_subpop(eigenvec, panel, study_ids, var_explained, output_dir)

    write_study_pca_table(eigenvec, study_ids, output_dir)

    print("\n--- ADMIXTURE ---")
    plot_admixture(results_dir, reference_dir, output_dir, panel, study_ids)

    print("\n" + "=" * 60)
    print(f"All plots saved to: {output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
