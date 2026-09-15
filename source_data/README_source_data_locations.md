# Source data — one folder per figure (mirrors 07_Code)

- **CD96 figures (Fig 5-7, S12-S26):** plot-level CSVs copied from the working package's `source_data/`,
  split by the figure that consumes them (traced from each script's `sd()` calls). These are the exact
  numeric values behind each plot; a Nature-style single Source Data file per figure can be assembled from them.
- **Legacy figures (Fig 1-4, S1-S11):** data files (<50 MB) copied from the original analysis tree,
  remapped to current figure numbers (see `../07_Code/README_code_locations.md` for the crosswalk).

## Files intentionally NOT copied (raw, repository-scale — deposit to GEO/SRA, not the Source Data files)
| File | Size | Figure | Legacy path |
|------|------|--------|-------------|
| `ALSF_RNA-seq_expr_HSPC_0.05_other_0.1.csv` | 3.7 GB | S10 screen | `ALSF_AML_plot/ALSF AML Figure7/` |
| `ALSF_RNA-seq_expr_HSPC_0.025/0.02/absent_other_0.1.csv` | 1.8-2.4 GB each | S10 screen | same |
| `ALSF_RNA-seq_Combo_obs.csv` | 63 MB | S10 | same |
| `receptor_ligand_interactions_mitab_v1.0_April2017.txt` | 84 MB | S10 (LR DB) | same |
| TARGET STAR counts / CPM / RPKM matrices | 143-576 MB | Fig 2/3/5 bulk | `ALSF_AML_plot/` root |
| `ALSF_AML_subtype_Count_matrix.tsv` | 576 MB | scRNA raw | `ALSF_AML_plot/` root |
| scRNA count matrix `count_matrix/ALSF_AML_raw_counts.h5ad`, SCENIC looms | GB | CD96 figs | working package / `ALSF_AML_plot/H5AD/` |

## Figures without a dedicated source-data table here
- **S1, S2** — rendered from the scanpy atlas; no separate numeric table (embeddings/marker dot plots).
- **S9, S11** — CITE-seq ADT data not located (gap; see MANIFEST).
