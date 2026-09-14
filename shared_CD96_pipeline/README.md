# shared_CD96_pipeline — self-contained figure bundle (CD96 story: Figs 5–7, deployed S14–S28)

Runnable from this folder with no external paths. Deployed supplement numbers are S14–S28; note
that the individual `figS_*.R` scripts still write output files on a superseded supplement
numbering (`RENAME_LOG.json` reconciles them to the deployed S-numbers). Main figure numbers match
the manuscript.

## Layout
- `scripts/`      — all R + Python code (`_style.R` defines colours and the `sd()` data loader);
  subfolders `cellchat/`, `upstream_builders/`, `legacy/`.
- `source_data/`  — the 35 plot-level CSVs the R renderers read (`sd("figNX_*.csv")`). Names match
  the manuscript figures; see `RENAME_LOG.json` for the historical old→new map.
- `figures/`      — output directory (created on run).
- `run_all.sh`    — master runner. `bash run_all.sh` renders Figs 5–7 + CD96 supplements from
  `source_data/`, assembles multi-panel figures, and (if bundled) rebuilds the combined PDF.

## Run
```bash
bash run_all.sh                      # uses Rscript + python3 on PATH
RSCRIPT=/path/to/Rscript PY=/path/to/python bash run_all.sh   # override interpreters
```
The R figure-rendering scripts need only the tidyverse/ggplot stack (see `REPRODUCIBILITY_CHECKLIST.md`
for exact packages + versions) and read only from `source_data/`. The `scripts/upstream_builders/`
(Seurat / CZ CELLxGENE Census / cohort matrices) regenerate `source_data/` from raw inputs and are
**not** needed to reproduce the figures; their large raw inputs are deposited to GEO.

## Provenance
- `RENAME_LOG.json` — old→new source-data filename map (reverse it to undo the renumbering).
- `SOURCE_DATA.md` — per-figure script → source-data → output map.
