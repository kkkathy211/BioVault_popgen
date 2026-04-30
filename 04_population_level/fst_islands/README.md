# fst_islands — Pairwise FST Between Caribbean Islands

Compute Wright's fixation index (FST) for every pair of Caribbean islands, then visualize the resulting matrix as a heatmap, hierarchical clustering dendrogram, and PCA.

## Inputs

```
data/ (Did not upload here!)
├── allele_freq_barbados.tsv     # per-locus AF for Barbados
├── allele_freq_bahamas.tsv      # Bahamas
├── allele_freq_BVI.tsv          # British Virgin Islands
├── allele_freq_bermuda.tsv      # Bermuda
├── allele_freq_stlucia.tsv      # St. Lucia
└── allele_freq_TT.tsv           # Trinidad & Tobago
```

Each file is tab-separated with columns: `locus_key`, `allele_frequency`, `allele_count`. ~40-50 MB per file. **Gitignored** — provided by collaborators (not regeneratable from public data within this repo).

## Scripts (run in order)

```
scripts/
├── 01_load_merge.py        # Read 6 island TSVs, intersect on shared loci, write merged AF table
├── 02_compute_fst.py       # Pairwise FST matrix (Hudson, Weir-Cockerham, or similar)
├── 03_visualize.py         # Heatmap, clustermap, dendrogram, PCA
├── figure3_publishable.py  # Publication-grade Figure 3
└── run_pipeline.sh         # Orchestrator
```

## Run

```bash
cd scripts
bash run_pipeline.sh
```

## Outputs

```
results/
├── merged/                # Intersected per-locus AF table (all 6 islands)
├── fst/                   # Pairwise FST matrix (TSV)
└── pca/                   # PCA on the AF matrix (eigenvec, eigenval)

plots/
├── fst_heatmap.png
├── fst_heatmap_clustermap.png   # row/col-clustered heatmap
├── fst_dendrogram.png           # hierarchical clustering tree
├── pca_populations.png          # populations in PC space
├── figure2_C_D.png              # composite figure
└── figure3_publishable.pdf      # publication-grade vector

logs/                            # one log per pipeline stage
```
