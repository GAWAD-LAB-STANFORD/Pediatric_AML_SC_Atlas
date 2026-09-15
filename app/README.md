# Pediatric AML Atlas Browser

A Shiny browser over the manuscript Source Data: surface-target ranking, per-leukemic-state
expression, normal-tissue safety, and combination coverage. It reads only the CSVs in
`../source_data`, so it runs without the atlas `.h5ad`.

```r
install.packages(c("shiny","bslib","ggplot2","DT","dplyr","hexbin"))
shiny::runApp("app")        # from the repository root
```

## Panels

| Panel | What it shows | Source |
| :-- | :-- | :-- |
| **Targets** | All 1,704 screened surfaceome genes ranked by the composite score, with % leukemic cells, patient breadth, and normal HSPC/myeloid expression | `Figure_5/fig5A_top15.csv` |
| **Leukemic states** | Percent positive across the leukemic states, coloured by each state's TARGET outcome association, with the normal-HSPC baseline | `Figure_5/fig5D_ls.csv` |
| **Safety** | Normal-marrow cell types and vital-organ expression for the twelve antigens carried past the Figure 5 screen | `Figure_6/fig6tox_hema.csv`, `fig6tox_organ.csv` |
| **Combinations** | Every single antigen and "X OR Y" pair, plotted as patient coverage against normal-marrow toxicity | `Figure_7/fig7_allpairs.csv` |
| **Atlas** | Single-cell UMAP, by cell density or Monocle pseudotime | `Figure_3/v3/monocle_umap_coords.csv` |

Selecting a gene — from the dropdown or by clicking a row in the Targets table — drives every panel.
Genes that were screened but not carried forward (most of the 1,704) show the ranking panels and a
note explaining that no safety or combination data exists for them.

## Notes

The atlas panel bins cells into hexagons, so rendering cost scales with the grid rather than the cell
count — the same approach used in the lab's earlier
[CellSeek](https://github.com/GAWAD-LAB-STANFORD/CellSeek) explorer. The UMAP file is subsampled to
1,000 cells per group for display.
