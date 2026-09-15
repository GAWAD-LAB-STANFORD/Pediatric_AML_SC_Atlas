# Pediatric AML Target Discovery

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
| **Discover** | Interactive screen: set your own efficacy, breadth and marrow-sparing thresholds, reweight the composite score, watch the funnel recompute, and export the passing genes | `Figure_5/fig5A_top15.csv`, `fig5A_funnel.csv` |
| **Target profile** | Per-leukemic-state expression coloured by outcome association, position against the 48 head-to-head antigens, and coverage across 14 cytogenetic subtypes | `Figure_5/fig5D_ls.csv`, `fig5C_headtohead.csv`, `fig5D_cyto.csv` |
| **Safety** | Normal-marrow cell types, vital-organ expression, and the five-axis profile for the twelve antigens carried past the screen | `Figure_6/fig6tox_hema.csv`, `fig6tox_organ.csv`, `fig6D_radar.csv` |
| **Validation** | Whether a target holds beyond discovery: paediatric TARGET bulk, adult Beat AML bulk, and an independent adult single-cell cohort | `Figure_5/fig5F_meta_all.csv`, `fig5F_adult_beataml_meta.csv`, `Figure_S23/adult_scrna_overall.csv` |
| **Combinations** | Every single antigen and "X OR Y" pair, plotted as patient coverage against normal-marrow toxicity | `Figure_7/fig7_allpairs.csv` |
| **Atlas** | Monocle trajectory embedding of the **paediatric** cohort (49 leukemic states + 8 normal marrow compartments), by compartment, pseudotime, or leukemic fraction | `Figure_3/v3/monocle_umap_coords.csv` |

Selecting a gene — from the dropdown or by clicking a row in the Targets table — drives every panel.
Genes that were screened but not carried forward (most of the 1,704) show the ranking panels and a
note explaining that no safety or combination data exists for them.

## Which cohort is which

Everything is the **paediatric** discovery cohort except the **Validation** panel, which is the
explicit adult comparison (Beat AML bulk and an independent adult single-cell cohort). The atlas
embedding is paediatric: its 57 groups are the 49 leukemic states `LS_1`-`LS_49` plus eight
`NORM_*` marrow compartments.

## Notes

Each of the 57 groups in the embedding is subsampled to exactly 1,000 cells, so **apparent cell
density reflects that subsampling rather than true abundance** — compartment and pseudotime are the
meaningful readouts, and a raw density view is deliberately not offered. Hex-binned panels follow the
approach used in the lab's earlier [CellSeek](https://github.com/GAWAD-LAB-STANFORD/CellSeek) explorer,
where render cost scales with the grid rather than the cell count.

Only a subset of the 1,704 screened genes carries the deeper panels (19 have per-state and subtype
data, 12 have safety, validation and combination data). Panels without data for the selected gene say
so explicitly rather than rendering an empty plot.
