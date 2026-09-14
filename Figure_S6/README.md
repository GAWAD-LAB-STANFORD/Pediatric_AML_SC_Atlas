# Figure S6 — Annotation of the 49 leukemic states

**Self-contained.** Reproduces deployed **Figure S7** (`../../05_Supplementary_Figures/Figure_S6__leukemic_state_annotation.pdf`) from its own source data.
Part of the atlas/prognosis analysis stream (legacy Figure_1-4 pipelines); scripts here are copies so this figure stands alone.

## Reproduce
```bash
cd scripts
Rscript panel_state_annotation.R
```
Reads `../source_data/` and writes `../figures/`.

## Contents
- `scripts/` — compute_state_annotation.py, panel_state_annotation.R + shared manuscript_palette.R/config.py
- `source_data/` — 2 file(s): state_annotation_heatmap_z.csv, state_annotation_meta.csv
