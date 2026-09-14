# Figure S3 — Selection of the 49 reproducible leukemic states

**Self-contained.** Reproduces deployed **Figure S3** (`../../05_Supplementary_Figures/Figure_S3__leukemic_state_selection.png`) from its own source data.
Part of the atlas/prognosis analysis stream (legacy Figure_1-4 pipelines); scripts here are copies so this figure stands alone.

## Reproduce
```bash
cd scripts
Rscript panel_state_selection.R
```
Reads `../source_data/` and writes `../figures/`.

## Contents
- `scripts/` — panel_state_selection.R + shared manuscript_palette.R/config.py
- `source_data/` — 1 file(s): jaccard_percluster.csv
