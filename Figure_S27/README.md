# Figure S27 — Every two-target combination landscape

**Self-contained.** This folder reproduces deployed **Figure S27** (`../../05_Supplementary_Figures/Figure_S27__allpairs_window.pdf`) from its own
plot-level source data. Part of the integrated CD96 analysis (canonical pipeline:
`../shared_CD96_pipeline/`, `run_all.sh`); the scripts here are copies so this figure stands alone.

## Reproduce (from `scripts/`, given the R/Python environment in `../shared_CD96_pipeline/REPRODUCIBILITY_CHECKLIST.md`)
```bash
cd scripts
Rscript fig7_allpairs.R
```
The renderer reads `../source_data/*.csv` and writes `figures/`. Deployed PDF is in `05_Supplementary_Figures/` (or `04_Main_Figures/`).

## Contents
- `scripts/` — renderer(s): fig7_allpairs.R; builder(s): build_fig7_allpairs.py; shared `_style.R` + `config.py` (bundled)
- `source_data/` — plot-level source data (1 file): fig7_allpairs.csv

## Notes
- Renderers reproduce the figure from `source_data/` alone (no raw data needed).
- Builders (`build_*.py`) regenerate `source_data/` from the raw single-cell atlas / TARGET bulk
  (GEO/dbGaP; path in `config.py` `H5`), so they require the deposited raw data to re-run.
