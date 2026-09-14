# Figure S9 — Transcription-factor regulon detail

Reproduces deployed **Figure S9** (`../../05_Supplementary_Figures/Figure_S9__regulon_detail.pdf`).
Part of the atlas/prognosis analysis stream (legacy Figure_1–4 pipelines); scripts here are copies so this figure stands alone.

## Panels
- **A** — per-cell-type regulon rank plots (SCENIC/AUCell).
- **B** — population × regulon activity heatmap (canonical healthy-BM regulons recovered as validation).
- **C** — full 57-state × regulon activity landscape (the homeobox/HOX subset of this landscape is Figure **2C**). Added 2026-08.

## Reproduce
```
../shared_CD96_pipeline/.venv/bin/python scripts/composite_S9.py
```
`composite_S9.py`:
1. renders the external pySCENIC A/B source (`./Figure_S9__regulon_detail_AB_source.pdf`, kept in this dir),
2. patches the stale baked-in "Figure S5" title → "Figure S9" (the file predates the S5→S9 supplementary renumber; the number is baked into the PDF, not added at layout),
3. appends panel **C** = `../../04_Main_Figures/Figure_2_panels/Figure_2_regulon_states.png` (built by `../Figure_2/panelA2_regulon_states.py`),
and writes the deployed `Figure_S9__regulon_detail.{png,pdf}`. Rerun after regenerating either input.

## Inputs
- `Figure_S9__regulon_detail_AB_source.pdf` (this dir) — external pySCENIC A/B render (title-patched at composite time; do not delete).
- Panel C's matrix: `../../08_Source_Data/Figure_S9/regulon_landscape_z.csv` (+ `LS_annotation.csv`, `regulon_landscape_mean.csv`, `proliferation_strip.csv`, `LS_wholecohort_prognosis.csv`) — copies of the canonical regulon-compute outputs in `../../08_Source_Data/Figure_3/v3/`; regenerate via `../Figure_3/compute_regulon_landscape.py`. This is the same landscape Figure 2C draws its HOX subset from.

## Source data
Canonical plot-level source data is in **`../../08_Source_Data/Figure_S9/`** (consistent with the other figures — the code tree holds no `source_data/` copy). The pySCENIC TF lists + top-regulon tables (`312TF_ALSF_AML.txt`, `topreg_*.csv`) live there too. A/B are rendered from external pySCENIC (AUCell) output; .loom outputs go to GEO (GAP: deposit).

## Contents
- `scripts/composite_S9.py` — in-tree renderer (A/B source + title patch + panel C).
- `Figure_S9__regulon_detail_AB_source.pdf` — external A/B render (composite input).
- `_superseded/` — archived redundant code-tree `source_data/` copy (data is canonical in `08_Source_Data/Figure_S9/`).
