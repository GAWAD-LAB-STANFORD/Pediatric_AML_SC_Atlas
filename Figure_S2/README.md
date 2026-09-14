# Figure S2 — AML atlas covariates

**Self-contained.** Reproduces deployed **Figure S2** (`../../05_Supplementary_Figures/Figure_S2__AML_atlas_covariates.png`) from its own source data.
Part of the atlas/prognosis analysis stream (legacy Figure_1-4 pipelines); scripts here are copies so this figure stands alone.

## Reproduce
```bash
cd scripts
python panel_S2_composite.py
```
Reads `../source_data/` and writes `../figures/`.

## Note
Renders covariate UMAPs directly from the scanpy atlas (config H5); needs the GEO atlas to re-run.

## Contents
- `scripts/` — panel_S2_composite.py
- `source_data/` — 1 file(s): README_Figure_S2.md
