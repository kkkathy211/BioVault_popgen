# 04 — Population-Level Analyses

Analyses whose primary input is **per-population allele frequency tables**, not per-individual genotypes. Each row is a SNP; each column (or each file) is a population.

## The two analyses

| Sub-analysis | What it does | Input |
|---|---|---|
| [`fst_islands/`](fst_islands/) | Compute pairwise FST between Caribbean islands and visualize as heatmap, dendrogram, and PCA. | Per-island allele frequency TSVs (Barbados, Bahamas, BVI, Bermuda, St. Lucia, Trinidad & Tobago) |
| [`aims_differential_snps/`](aims_differential_snps/) | Identify ancestry-informative markers (AIMs) and per-island differential SNPs by comparing each island's allele frequencies against gnomAD reference populations. | Per-island AF + gnomAD v3 sites AF |

## Independence

These two analyses are independent. They take different input formats (one needs the island-level TSVs from collaborators; the other needs gnomAD downloads) and produce different output types.
