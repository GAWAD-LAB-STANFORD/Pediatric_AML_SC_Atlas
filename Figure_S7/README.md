# Figure S7 — State-level stemness does not distinguish prognosis

**Self-contained.** Reproduces deployed **Figure S9** (`../../05_Supplementary_Figures/Figure_S7__stemness_by_prognosis.png`) from its own source data.
Part of the atlas/prognosis analysis stream (legacy Figure_1-4 pipelines); scripts here are copies so this figure stands alone.

## Reproduce
```bash
cd scripts
Rscript panelD_stem_by_prog.R
```
Reads `../source_data/` and writes `../figures/`.

## Contents
- `scripts/` — panelD_build_stem.py, panelD_stem_by_prog.R + shared manuscript_palette.R/config.py
- `source_data/` — 2 file(s): LS_stem_percell.csv, LS_wholecohort_prognosis.csv
