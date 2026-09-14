# Figure 5 — CD96 target discovery

**Self-contained.** This folder reproduces deployed **Figure 5** (`../04_Main_Figures/Figure_5__CD96_target_discovery.pdf`) from its own
plot-level source data. Part of the integrated CD96 analysis (canonical pipeline:
`../shared_CD96_pipeline/`, `run_all.sh`); the scripts here are copies so this figure stands alone.

## Reproduce (from `scripts/`, given the R/Python environment in `../shared_CD96_pipeline/REPRODUCIBILITY_CHECKLIST.md`)
```bash
cd scripts
Rscript fig5.R   # panels A-F
```
The renderer reads `../source_data/*.csv` and writes `figures/`. Deployed PDF is in `05_Supplementary_Figures/` (or `04_Main_Figures/`).

## Contents
- `scripts/` — renderer(s): fig5.R; builder(s): build_adult_bulk.py, build_fig5A_funnel.py, build_fig5B_umaps.py, build_fig5CD_leads.py, build_fig5C_headtohead.py, build_fig5D_all.py, build_fig5D_ls.py, build_fig5D_panel.py, build_fig5E_target_leads.py; shared `_style.R` + `config.py` (bundled)
- `source_data/` — plot-level source data (9 files): fig5A_funnel.csv, fig5A_top15.csv, fig5C_headtohead.csv, fig5D_cyto.csv, fig5D_ls.csv, fig5F_adult_beataml_meta.csv, fig5F_adult_beataml_violin.csv, fig5F_meta_all.csv, fig5F_violin.csv

## Notes
- Renderers reproduce the figure from `source_data/` alone (no raw data needed).
- Builders (`build_*.py`) regenerate `source_data/` from the raw single-cell atlas / TARGET bulk
  (GEO/dbGaP; path in `config.py` `H5`), so they require the deposited raw data to re-run.
