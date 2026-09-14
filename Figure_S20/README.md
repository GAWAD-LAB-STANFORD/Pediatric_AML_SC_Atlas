# Figure S20 — Candidate antigens by cytogenetic subtype, pediatric TARGET

**Self-contained.** This folder reproduces deployed **Figure S21** (`../../05_Supplementary_Figures/Figure_S20__cyto_violins_pediatric.pdf`) from its own
plot-level source data. Part of the integrated CD96 analysis (canonical pipeline:
`../shared_CD96_pipeline/`, `run_all.sh`); the scripts here are copies so this figure stands alone.

## Reproduce (from `scripts/`, given the R/Python environment in `../shared_CD96_pipeline/REPRODUCIBILITY_CHECKLIST.md`)
```bash
cd scripts
Rscript figS_cyto_violin.R
```
The renderer reads `../source_data/*.csv` and writes `figures/`. Deployed PDF is in `05_Supplementary_Figures/` (or `04_Main_Figures/`).

## Contents
- `scripts/` — renderer(s): figS_cyto_violin.R; builder(s): build_figS_cyto_pediatric.py; shared `_style.R` + `config.py` (bundled)
- `source_data/` — plot-level source data (3 files): fig5F_meta_all.csv, figS_cyto_adult_violin.csv, figS_cyto_pediatric_violin.csv

## Notes
- Renderers reproduce the figure from `source_data/` alone (no raw data needed).
- Builders (`build_*.py`) regenerate `source_data/` from the raw single-cell atlas / TARGET bulk
  (GEO/dbGaP; path in `config.py` `H5`), so they require the deposited raw data to re-run.
