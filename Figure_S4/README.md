# Figure S4 — Bulk-deconvolution validation

**Self-contained.** Reproduces deployed **Figure S4** (`../../05_Supplementary_Figures/Figure_S4__bulk_deconvolution.pdf`) from its own source data.
Part of the atlas/prognosis analysis stream (legacy Figure_1-4 pipelines); scripts here are copies so this figure stands alone.

## Reproduce
```bash
cd scripts
Rscript panel_deconv_validation.R
```
Reads `../source_data/` and writes `../figures/`.

## Contents
- `scripts/` — compute_deconv_qc.py, panel_deconv_validation.R + shared manuscript_palette.R/config.py
- `source_data/` — 3 file(s): LS_detection_49.csv, LS_wholecohort_prognosis.csv, deconv_fit_qc.csv
