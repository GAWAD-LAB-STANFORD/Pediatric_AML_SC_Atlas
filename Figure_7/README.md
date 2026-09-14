# Figure 7 — Combination strategy

**Self-contained.** This folder reproduces deployed **Figure 7** (`../04_Main_Figures/Figure_7__combination_strategy.pdf`) from its own
plot-level source data. Part of the integrated CD96 analysis (canonical pipeline:
`../shared_CD96_pipeline/`, `run_all.sh`); the scripts here are copies so this figure stands alone.

## Reproduce (from `scripts/`, given the R/Python environment in `../shared_CD96_pipeline/REPRODUCIBILITY_CHECKLIST.md`)
```bash
cd scripts
Rscript cellchat/fig7_circos_3panel.R && Rscript fig7_combo.R && Rscript fig7d_bestcombo.R && Rscript fig7_allpairs.R && python build_fig7_assemble.py
```
The renderer reads `../source_data/*.csv` and writes `figures/`. Deployed PDF is in `05_Supplementary_Figures/` (or `04_Main_Figures/`).

## Contents
- `scripts/` — renderer(s): cellchat_CD96.R, fig7_allpairs.R, fig7_circos_3panel.R, fig7_combo.R, fig7d_bestcombo.R; builder(s): build_cellchat_3path_input.py, build_fig7_allpairs.py, build_fig7_assemble.py, build_fig7_combo.py; shared `_style.R` + `config.py` (bundled)
- `source_data/` — plot-level source data (2 files): fig7_allpairs.csv, fig7_combo.csv

## Notes
- Renderers reproduce the figure from `source_data/` alone (no raw data needed).
- Builders (`build_*.py`) regenerate `source_data/` from the raw single-cell atlas / TARGET bulk
  (GEO/dbGaP; path in `config.py` `H5`), so they require the deposited raw data to re-run.
