# Figure S5 — Whole-cohort EFS hazard ratios for the 49 states

**Self-contained.** Reproduces deployed **Figure S8** (`../../05_Supplementary_Figures/Figure_S5__state_EFS_forest.png`) from its own source data.
Part of the atlas/prognosis analysis stream (legacy Figure_1-4 pipelines); scripts here are copies so this figure stands alone.

## Reproduce
```bash
cd scripts
Rscript panel_state_hr_forest.R
```
Reads `../source_data/` and writes `../figures/`.

## Contents
- `scripts/` — km_table_all_states.R, panel_state_hr_forest.R + shared manuscript_palette.R/config.py
- `source_data/` — 2 file(s): LS_annotation.csv, LS_wholecohort_cox_efs.csv
