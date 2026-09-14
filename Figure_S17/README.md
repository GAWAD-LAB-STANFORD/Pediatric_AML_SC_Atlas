# Figure S17 — CD96-comparator co-expression (metacells + bulk)

**Self-contained.** This folder reproduces deployed **Figure S17** (`../../05_Supplementary_Figures/Figure_S17__correlations.pdf`) from its own
plot-level source data. Part of the integrated CD96 analysis (canonical pipeline:
`../shared_CD96_pipeline/`, `run_all.sh`); the scripts here are copies so this figure stands alone.

## Reproduce (from `scripts/`, given the R/Python environment in `../shared_CD96_pipeline/REPRODUCIBILITY_CHECKLIST.md`)
```bash
cd scripts
Rscript figS.R   # correlations section
```
The renderer reads `../source_data/*.csv` and writes `figures/`. Deployed PDF is in `05_Supplementary_Figures/` (or `04_Main_Figures/`).

## Contents
- `scripts/` — renderer(s): figS.R; builder(s): build_FigureS_CD96_correlations.py, build_bulk_S17.py, build_metacells_S7.py; shared `_style.R` + `config.py` (bundled)
- `source_data/` — plot-level source data (6 files): fig5F_meta_all.csv, fig5F_violin.csv, figS17_corr_r.csv, figS17_points_bulk.csv, figS17_points_singlecell.csv, figS_cyto_pediatric_violin.csv

## Notes
- Renderers reproduce the figure from `source_data/` alone (no raw data needed).
- Builders (`build_*.py`) regenerate `source_data/` from the raw single-cell atlas / TARGET bulk
  (GEO/dbGaP; path in `config.py` `H5`), so they require the deposited raw data to re-run.
