# Figure S12 — CITE-seq surface-protein validation of the T/NK subsets

**Self-contained.** Reproduces deployed **Figure S12** (`../../05_Supplementary_Figures/Figure_S12__T_NK_CITEseq_surface.pdf`) from its own source data.
Part of the atlas/prognosis analysis stream (legacy Figure_1-4 pipelines); scripts here are copies so this figure stands alone.

## Reproduce
```bash
cd scripts
python build_figS11_citeseq.py
```
Reads `../source_data/` and writes `../figures/`.

## Note
Builder reads the CITE-seq object (config; GEO) and maps ADT onto the T/NK embedding.

## Contents
- `scripts/` — build_figS11_citeseq.py
- `source_data/` — 3 file(s): README_Figure_S11.md, figS11_surface_pct_by_subset.csv, figS11_tnk_surface_percell.csv
