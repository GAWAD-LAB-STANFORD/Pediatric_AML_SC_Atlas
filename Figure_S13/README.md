# Figure S13 — T-cell receptor repertoire of the marrow T/NK subsets

**Self-contained.** Reproduces deployed **Figure S13** (`../../05_Supplementary_Figures/Figure_S13__TNK_TCR_repertoire.pdf`) from its own source data.
Part of the atlas/prognosis analysis stream (legacy Figure_1-4 pipelines); scripts here are copies so this figure stands alone.

## Reproduce
```bash
cd scripts
Rscript panel_tnk_tcr_extra.R && python composite_supp_tcr.py
```
Reads `../source_data/` and writes `../figures/`.

## Note
Renderer panels build from the shipped TCR source data; compute_tnk_tcr2.py regenerates them from the GEO VDJ object. composite stitches the rendered panels.

## Contents
- `scripts/` — composite_supp_tcr.py, compute_tnk_tcr2.py, panel_tnk_tcr_extra.R + shared manuscript_palette.R/config.py
- `source_data/` — 6 file(s): tcr_clonality_by_cyto.csv, tcr_clonality_downsampled.csv, tcr_exhaustion_by_expansion.csv, tcr_sharing.csv, tcr_trbv_usage.csv, tnk_subset_survival.csv
