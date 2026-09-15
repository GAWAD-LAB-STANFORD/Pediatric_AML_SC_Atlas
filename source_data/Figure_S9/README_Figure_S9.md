# Figure S9 — Transcription-factor regulon detail (source data)

Plot-level source data for deployed **Figure S9** (renderer: `../../07_Code/Figure_S9/scripts/composite_S9.py`).

## Panels A/B — per-cell-type regulon rank plots + population heatmap (external pySCENIC)
- `312TF_ALSF_AML.txt` — 312-TF input list for the pySCENIC run.
- `topreg_enriched Leukemia_with_normal_BM.csv`, `topreg_normal_BM.csv` — top-regulon summary tables.
- Rank plots/heatmap come from the external pySCENIC run (see `../Figure_3/`). GAP: deposit the full pySCENIC .loom output (to GEO).

## Panel C — full 57-state × regulon activity landscape (added 2026-08)
Plot-level matrices behind panel C (built by `../../07_Code/Figure_2/panelA2_regulon_states.py`). **Copies of the canonical regulon-compute outputs in `../Figure_3/v3/`** (shipped here so the figure is self-contained; regenerate from `../../07_Code/Figure_3/compute_regulon_landscape.py`):
- `regulon_landscape_z.csv` — z-scored AUCell regulon activity, 57 columns (49 leukemic states + 8 normal compartments) × regulons. Same matrix underlies the homeobox/HOX subset in **Figure 2C**.
- `regulon_landscape_mean.csv` — mean AUCell (for dot sizing).
- `proliferation_strip.csv` — per-column proliferation-module z (top annotation).
- `LS_annotation.csv` — per-state lineage-program labels.
- `LS_wholecohort_prognosis.csv` — per-state EFS prognosis (column colours); canonical copy in `../Figure_2/v3/`.
