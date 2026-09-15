# Code — one self-contained folder per figure

`07_Code/` has one folder per deployed figure (main Figures 1–7, Supplementary Figures **S1–S27**),
and `source_data/` mirrors the same numbering. Folder numbers now match the deployed figures
(the earlier superseded/drifted supplement numbering has been reconciled).

## Self-contained per-figure folders
Each CD96 figure (Figs 5–7) and every supplement folder (`Figure_S1/…Figure_S27/`) contains its own
runnable code and data, so a reviewer can reproduce that figure from its own folder:
```
Figure_<N>/
  scripts/       renderer(s) + builder(s) + bundled shared deps (_style.R / config.py / manuscript_palette.R)
  source_data/   the plot-level CSVs the renderer reads
  figures/       render output (the deployed copy is in 04_Main_Figures / 05_Supplementary_Figures)
  README.md      one-figure reproduce guide
```
Reproduce a figure: `cd Figure_<N>/scripts && Rscript <renderer>.R` (or `python <renderer>.py`), which
reads `../source_data/` and writes `../figures/`. Builders (`build_*.py`) regenerate `source_data/`
from the raw atlas / TARGET bulk (GEO/dbGaP; path in `config.py`) and need the deposited raw data.

Legacy main figures `Figure_1/…Figure_4/` keep their original multi-script pipelines (they already
ship the full scripts; their plot-level data is in `source_data/Figure_1…4/`).

## Two provenance streams
- **Atlas / prognosis (Figs 1–4, Suppl. S1–S13)** — leukemic-state atlas, prognosis, T/NK and CITE-seq
  validation. Built by the per-figure folders; three figures come from external pipelines
  (S1 scanpy, S5 pySCENIC, S11 muon/ADT) whose plot-level source data is shipped and whose raw
  objects deposit to GEO.
- **CD96 target discovery (Figs 5–7, Suppl. S14–S27)** — the 2026 revision. The canonical integrated
  pipeline is kept whole in **`shared_CD96_pipeline/`** (`bash run_all.sh` rebuilds every CD96 figure
  from bundled source data); the per-figure CD96 folders are self-contained copies of its scripts + data.

## Current ← original main-figure crosswalk
| Current | Original | Content |
|---|---|---|
| Fig 1–3 | Fig 1–3 | atlas / 49 leukemic states + prognosis / LS deep-dive |
| Fig 4 | **Fig 5** | T/NK-cell states + prognosis |
| Fig 5–7 | (new) | CD96 discovery / toxicity / combination |

## Supplement provenance (deployed S1–S27)
See `FIGURE_SCRIPT_MAP.md` for the full per-figure table (folder, renderer, provenance). Key points:
- **Fig-4 T/NK trio** is ordered normalization (**S11**) → validation (**S12**) → TCR repertoire (**S13**).
- **S27** is the every-two-target combination landscape (was S28; the numbering gap is closed).
- **Dropped/merged** (not in the deployed set, archived under `Archive_superseded/`): the surfaceome
  MCC/F1 screen (old S12) and the per-PPAC additive window (old S27, merged into S25); the two hdWGCNA
  figures (module–outcome result did not reproduce on the corrected 49-state scheme).

## Provenance gaps for deposit
- **SCENIC (S5)** — plot-level AUCell Z-scores + TF lists shipped; pySCENIC `.loom` outputs to GEO.
- **CITE-seq (S11 normalization / S12 validation)** — S12 is fully in-tree (`build_figS11_citeseq.py`);
  S11 ADT normalization ran in an external muon notebook; the CITE-seq/muon object deposits to GEO.
- **Atlas UMAPs (S1, and Fig 1 UMAPs)** — external scanpy pipeline; plot-level coordinates shipped.
- **TCR (S13)** — renderer reproduces from shipped TCR source data; the Cell Ranger VDJ object deposits
  to GEO for a full rebuild via `compute_tnk_tcr2.py`.

See `../MANIFEST.md` for the full inventory.
