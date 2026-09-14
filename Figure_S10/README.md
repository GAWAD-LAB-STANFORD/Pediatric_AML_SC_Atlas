# Figure S10 — Cellular composition and cell cycle

**Self-contained.** Reproduces deployed **Figure S6** (`../../05_Supplementary_Figures/Figure_S10__composition_cell_cycle.png`) from its own source data.
Part of the atlas/prognosis analysis stream (legacy Figure_1-4 pipelines); scripts here are copies so this figure stands alone.

## Reproduce
```bash
cd scripts
Rscript panel_S9_composite.R
```
Reads `../source_data/` and writes `../figures/`.

## Contents
- `scripts/` — compute_S9_composition.py, panel_S9_composite.R + shared manuscript_palette.R/config.py
- `source_data/` — 4 file(s): s9_lineage_by_sample.csv, s9_normal_lineage_by_sample.csv, s9_phase_by_sample.csv, state_cellcycle_by_prognosis.csv
