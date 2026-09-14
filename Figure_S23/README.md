# Figure S23 — Adult AML single-cell target validation

**Self-contained.** This folder reproduces deployed **Figure S23** (`../../05_Supplementary_Figures/Figure_S23__adult_scrna.pdf`) from its own
plot-level source data. Part of the integrated CD96 analysis (canonical pipeline:
`../shared_CD96_pipeline/`, `run_all.sh`); the scripts here are copies so this figure stands alone.

## Reproduce (from `scripts/`, given the R/Python environment in `../shared_CD96_pipeline/REPRODUCIBILITY_CHECKLIST.md`)
```bash
cd scripts
Rscript fig_adult_scrna.R
```
The renderer reads `../source_data/*.csv` and writes `figures/`. Deployed PDF is in `05_Supplementary_Figures/` (or `04_Main_Figures/`).

## Contents
- `scripts/` — renderer(s): fig_adult_scrna.R; builder(s): build_adult_scrna_markers.py, build_adult_scrna_umap.py; shared `_style.R` + `config.py` (bundled)
- `source_data/` — plot-level source data (2 files): adult_scrna_overall.csv, adult_scrna_subtype.csv

## Notes
- Renderers reproduce the figure from `source_data/` alone (no raw data needed).
- Builders (`build_*.py`) regenerate `source_data/` from the raw single-cell atlas / TARGET bulk
  (GEO/dbGaP; path in `config.py` `H5`), so they require the deposited raw data to re-run.
