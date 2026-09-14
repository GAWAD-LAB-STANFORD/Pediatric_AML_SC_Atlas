# Figure 6 — Normal-tissue toxicity

**Self-contained.** This folder reproduces deployed **Figure 6** (`../04_Main_Figures/Figure_6__normal_tissue_toxicity.pdf`) from its own
plot-level source data. Part of the integrated CD96 analysis (canonical pipeline:
`../shared_CD96_pipeline/`, `run_all.sh`); the scripts here are copies so this figure stands alone.

## Reproduce (from `scripts/`, given the R/Python environment in `../shared_CD96_pipeline/REPRODUCIBILITY_CHECKLIST.md`)
```bash
cd scripts
Rscript fig6_toxicity.R && Rscript fig6_DE.R && python build_fig6_assemble.py
```
The renderer reads `../source_data/*.csv` and writes `figures/`. Deployed PDF is in `05_Supplementary_Figures/` (or `04_Main_Figures/`).

## Contents
- `scripts/` — renderer(s): fig6_DE.R, fig6_toxicity.R; builder(s): build_fig6A_spec_leads.py, build_fig6_DE.py, build_fig6_assemble.py, build_fig6_toxicity.py, query_census_normal_bodywide.py; shared `_style.R` + `config.py` (bundled)
- `source_data/` — plot-level source data (4 files): fig6D_radar.csv, fig6tox_celltypes.csv, fig6tox_hema.csv, fig6tox_organ.csv

## Notes
- Renderers reproduce the figure from `source_data/` alone (no raw data needed).
- Builders (`build_*.py`) regenerate `source_data/` from the raw single-cell atlas / TARGET bulk
  (GEO/dbGaP; path in `config.py` `H5`), so they require the deposited raw data to re-run.
