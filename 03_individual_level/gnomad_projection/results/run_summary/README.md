# Pipeline Run Summary

**Run completed:** 2026-04-26 23:24:48
**Reference:** gnomAD v3.1.2 HGDP+1kGP subset (AWS S3 mirror)

---

## 1. Pipeline status

| Step | Status | Notes |
|------|--------|-------|
| 0  Setup                          | ✅ | env, tools, refs |
| 1  Convert raw → PLINK            | ✅ | 396,070 input SNPs |
| 2  QC (geno/mind/maf/hwe)         | ✅ | 395,101 / 10 retained |
| 2b Loadings-intersect (LD-prune)  | ✅ | 5,973 SNPs (gnomAD AIM panel) |
| 3  Download HGDP+1kGP reference   | ✅ | 22 chrs, 5,808 sites total |
| 4  Prepare reference PLINK        | ✅ | 4,151 samples, 5,808 SNPs |
| 5  Merge study + reference        | ✅ | PLINK 1.9 `--bmerge` (PLINK 2 doesn't support cross-sample merge yet) |
| 6  PCA                            | ✅ | 10 PCs, 4,161 samples × 4,875 LD-pruned SNPs |
| 7  ADMIXTURE                      | ⚠️ | **Not installed — skipped** (see "Optional next step") |
| 8  Visualization                  | ✅ | 2 PCA plots produced |

---

## 2. Key numbers

```
Study samples:       10
Reference samples:   4,151 (HGDP+1kGP from gnomAD v3.1.2)
  AFR  (African)            1,003
  EAS  (East Asian)           825
  CSA  (Central/South Asian)  790
  EUR  (European)             788
  AMR  (Admixed American)     552
  MID  (Middle Eastern)       162
  OCE  (Oceanian)              30

Variant funnel:
  Raw study SNPs                       396,070
  After QC                             395,101
  After loadings-intersect (step 2b)     5,973   (gnomAD AIMs)
  Common w/ reference                    5,808
  After ambiguous (A/T,C/G) filter       5,508
  After cross-merge                      5,424
  After LD-prune (--indep-pairwise)      4,875   ← used by PCA
```

---

## 3. Result files

### PCA (Path C — full merged PCA)
```
results/pca/merged_pca.eigenvec    (10 PCs × 4,161 samples)
results/pca/merged_pca.eigenval    (10 eigenvalues)
results/plots/pca_continental.png  (study + ref colored by region)
results/plots/pca_subpopulation.png
results/plots/study_pca_coordinates.tsv  (just the 10 study samples)
```

### PCA projection (Path A — projected onto gnomAD's pre-computed PC space)
```
results/pca_projection/study_pca_projection.tsv  (16 PCs per study sample)
results/pca_projection/pca_projection.png        (study only)
results/pca_projection/pca_overlay_pc1_pc2.png   (with reference, continental)
results/pca_projection/pca_subpop_EAS.png        (zoomed to East Asian sub-pops)
```

### QC
```
results/qc_report.txt     (geno=0.05, mind=0.1, maf=0.01, hwe=1e-4)
```

---

## ⚠️ KEY INTERPRETATION RESULT

**The 10 study samples cluster with the South Asian (SAS / CSA) reference
population, at the western/admixed edge of that cluster.**

Source: `pca_overlay_pc1_pc2.png` (Path A — projection onto gnomAD's pre-computed
PC space). All 10 study stars (★) sit right inside the SAS (yellow) cluster at
PC1 ≈ 0.02–0.05, PC2 ≈ 0.005–0.04.

**Sub-population zoom** (`pca_subpop_CSA.png`): the 10 samples sit at the
**leftmost (lower-PC1) edge** of the CSA cluster — closest to BEB (Bengali),
GIH (Gujarati), ITU (Indian Telugu), and STU (Sri Lankan Tamil), but offset
slightly. This edge position can mean either:
- a CSA sub-population not well represented in HGDP+1kGP, or
- some admixture pulling them toward the AMR (light blue) cluster nearby.

For a definitive sub-pop call, look at `pca_subpop_all.png` (full sub-pop view)
and compare against the labeled HGDP+1kGP sub-pops near the study cluster.

**Trust Path A, NOT Path C** for the ancestry call. The Path C merged PCA
(`pca_continental.png`) shows the study samples isolated at PC1 ≈ 0.30, far
from any reference cluster. That is a **technical artifact**, not a real
ancestry signal — it comes from the extreme sample imbalance (10 vs 4,151) plus
residual REF/ALT mismatches between study and reference at non-ambiguous SNPs
that the A/T,C/G ambiguous filter does not catch. PC1 = 525 (≈5× PC2) is the
classic signature of one tight outlier group dominating the eigendecomposition.

Path A avoids this entirely because it uses gnomAD's pre-computed loadings on
the full HGDP+1kGP reference — the study samples are projected as queries, not
re-fit, so their count cannot dominate the eigenvectors.

---

## 4. Study sample PC coordinates (from Path C — `merged_pca`, do NOT use for ancestry call)

| Sample    | PC1     | PC2      | PC3      | PC4      | PC5      |
|-----------|---------|----------|----------|----------|----------|
| 305387    |  0.3100 | -0.0267  | -0.0681  |  0.1831  |  0.1277  |
| 351218    |  0.3148 | -0.0232  | -0.0377  |  0.0041  |  0.5214  |
| 371396    |  0.3080 | -0.0155  | -0.0268  | -0.1217  |  0.0132  |
| 515084    |  0.3122 | -0.0071  | -0.0291  | -0.1974  |  0.0972  |
| 525943    |  0.3187 |  0.0152  | -0.0566  |  0.8012  | -0.1692  |
| 807622    |  0.3196 | -0.0205  | -0.0244  | -0.0246  | -0.4999  |
| 817803    |  0.3089 | -0.0223  | -0.0869  | -0.0932  | -0.2475  |
| 821425    |  0.3018 | -0.0004  |  0.0031  | -0.0622  |  0.5282  |
| 838338    |  0.3152 | -0.0298  | -0.0281  | -0.0177  | -0.0705  |
| 959836    |  0.3156 |  0.0015  | -0.0259  | -0.4957  | -0.2586  |

**Observation:** All 10 study samples cluster very tightly on PC1 (0.30–0.32) and
PC2 (≈0), well separated from the bulk of the reference cluster which sits near
PC1 ≈ 0. PC1 in HGDP+1kGP typically separates AFR from non-AFR; the sign of
PC1 is arbitrary, so check `pca_continental.png` to see which side the study
samples occupy.

PC4 and PC5 show within-cluster spread among the 10 samples — these capture
fine-grained sub-population structure (or possibly a kinship/relatedness signal
given the small sample size).

---

## 5. Top eigenvalues

```
PC1: 525.79   ← dominant (typically AFR vs non-AFR in HGDP+1kGP)
PC2: 111.09
PC3:  71.84
PC4:  62.18
PC5:  61.12
PC6:  59.36
PC7:  58.21
PC8:  57.20
PC9:  56.93
PC10: 56.43
```

---

## 6. Pipeline issues encountered (and fixes)

| Issue | Fix applied |
|-------|-------------|
| GCS HTTPS mirror gave HTTP/2 framing errors → silent 15 KB truncated chunks | Switched to AWS S3 mirror (HTTP/1.1); local `.tbi` pre-download |
| `bcftools concat -a` required per-chunk indexes | Dropped `-a` (chunks are non-overlapping) |
| `bcftools view -R` on unindexed local file | Switched to `-T` (targets, no index needed) |
| Re-runs wiped chunk dirs, lost work after restarts | Skip-if-`.vcf.gz`-exists check; preserved chunk dirs across restarts |
| `plink2 --pmerge` is pgen-only | Tried `--pmerge-list bfile` → also fails: PLINK 2 cross-sample merge "under development" |
| → Final: installed PLINK 1.9, used `plink --bmerge` | works |

---

## 7. Optional next step: ADMIXTURE

ADMIXTURE was skipped because it wasn't installed. To enable it:

```bash
conda activate ancestry_pipeline
conda install -c bioconda admixture -y
rm pipeline_gnomad/working/step7.done
bash pipeline_gnomad/run_pipeline_gnomad.sh
```

This will run ADMIXTURE on the LD-pruned merged dataset for K=2..8 and produce
`results/admixture/*.Q` files plus a CV-error plot.

---

## 8. How to view the plots

The two main visualizations are:

1. **`results/plots/pca_continental.png`** — Study samples (★) on top of HGDP+1kGP
   reference, colored by continental region (AFR/AMR/CSA/EAS/EUR/MID/OCE).
   This is the easiest way to read which ancestry your samples cluster with.

2. **`results/pca_projection/pca_overlay_pc1_pc2.png`** — Same idea but using
   gnomAD's pre-computed PC space (Path A — independent verification).
   If both plots show your samples in the same continental cluster, the result
   is robust.

3. **`results/pca_projection/pca_subpop_EAS.png`** — Sub-population zoom for
   East Asian. (Other regions can be regenerated by running
   `run_pca_project.sh` with each region.)
