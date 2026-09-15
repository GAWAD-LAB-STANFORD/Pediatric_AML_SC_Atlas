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
| **Safety** | Browse the **top 500** targets on the marrow-safety axes (% normal HSPC vs % myeloid progenitors, available for all 1,704), then the selected gene's marrow cell types, vital organs and five-axis profile | `Figure_5/fig5A_top15.csv`, `Figure_6/fig6tox_hema.csv`, `fig6tox_organ.csv`, `fig6D_radar.csv` |
| **Validation** | Whether a target holds beyond discovery: paediatric TARGET bulk, adult Beat AML bulk, and an independent adult single-cell cohort | `Figure_5/fig5F_meta_all.csv`, `fig5F_adult_beataml_meta.csv`, `Figure_S23/adult_scrna_overall.csv` |
| **Combinations** | Every single antigen and "X OR Y" pair, plotted as patient coverage against normal-marrow toxicity | `Figure_7/fig7_allpairs.csv` |
| **Atlas** | The Figure 2A UMAP of the full **paediatric** atlas — 96,627 cells (70,108 leukemic across 49 states, 26,519 normal marrow) — coloured by compartment, cell type, leukemic state, or density | `Figure_2/fig2A_umap_cells.csv` |

Selecting a gene — from the dropdown or by clicking a row in the Targets table — drives every panel.
Genes that were screened but not carried forward (most of the 1,704) show the ranking panels and a
note explaining that no safety or combination data exists for them.

## Which cohort is which

Everything is the **paediatric** discovery cohort except the **Validation** panel, which is the
explicit adult comparison (Beat AML bulk and an independent adult single-cell cohort). The Atlas tab
shows the Figure 2A UMAP of all 96,627 paediatric cells; 70,108 carry a leukemic-state assignment
(`LS_1`-`LS_49`) and the remaining 26,519 are normal marrow.

## Notes

`Figure_2/fig2A_umap_cells.csv` is built by joining the atlas UMAP coordinates to
`Figure_3/v3/reference_assignments.csv`, mapping `res4.50` clusters to `LS_1`-`LS_49` via the stable
clusters in `jaccard_percluster.csv` (mean Jaccard >= 0.35, ordered by size) — the same rule the
figure scripts use. Because every cell is present rather than a per-group subsample, the density view
is meaningful. The hex-binned view follows the approach of the lab's earlier
[CellSeek](https://github.com/GAWAD-LAB-STANFORD/CellSeek) explorer, where render cost scales with the
grid rather than the cell count.

Only a subset of the 1,704 screened genes carries the deeper panels (19 have per-state and subtype
data, 12 have vital-organ, validation and combination data). Panels without data for the selected gene
say so explicitly rather than rendering an empty plot.

The distinction matters most in the Safety tab. **Marrow** safety (% normal HSPC, % myeloid
progenitors) is measured for all 1,704 genes, so the top-500 browser is real data throughout.
**Vital-organ** profiling against the CZ CELLxGENE Census exists only for the twelve antigens carried
past the screen — the browser marks which those are rather than leaving it implied.
