# Figure 2 — Leukemic states, prognosis, and function (corrected 49-state scheme)

Built on the **corrected** leukemic set: Cell Type ∈ {AML, AML-MKI67, AML-PCNA, AML-CD1C}
(70,108 cells; the earlier build wrongly included ~3,500 normal cells via `AML Sub-Clusters`).
Re-clustered at Leiden res 4.5, per-cluster Jaccard ≥ 0.35 → **49 stable states LS_1..LS_49**
(sorted by size). TARGET deconvolution = CIBERSORTx Jobs 24–27 (1,787 patients). All prognosis
is **event-free survival (EFS)**. Whole cohort first (A–E), then core-binding-factor (F–H).

Shared: `manuscript_palette.R` (ltc palette). Data: `../../08_Source_Data/Figure_2/v3/` and
`../../08_Source_Data/Figure_3/v3/` (LS_survival, LS_annotation, meancpm, reference_assignments,
LS_upgenes_wilcoxon). Panels render to `../../04_Main_Figures/Figure_2_panels/`.

## Panel → render script → data
| panel | content | render | data / compute |
|---|---|---|---|
| A | CD34 + 49 leukemic states on the atlas UMAP | `panelA_states_umap.py` | reference_assignments + h5ad |
| B | states by whole-cohort EFS prognosis (+ marker-gated normal-HSC reference) | `panelB_prognosis_umap.py` | `LS_wholecohort_prognosis.csv`, `normal_HSC_gated_barcodes.csv` |
| C | prognostic states × cytogenetic \| molecular + stemness strip | `panelC_states_subtypes.R` | LS_survival + molecular table; enrichment computed in-script |
| D | state-level stemness (LSC-Eppert, HSC) by prognosis group | `panelD_build_stem.py` → `panelD_stem_by_prog.R` | precomputed h5ad `LSC_EPPERT_Score`/`HSC_Score` → `LS_stem_percell.csv` |
| E | GO:BP of per-state up-genes (equal-N) + stemness strip | `panelE_prognostic_GO.R` | `LS_upgenes_wilcoxon.csv` (per-state Wilcoxon DE) |
| F | where CBF maps (t(8;21), inv(16)) on the UMAP | `panelF_cbf_umap.py` | reference_assignments (Cytogenetic) |
| G | within-CBF EFS-prognostic states (Cox forest) | `panelF_cbf_states.R` | LS_survival (CBF subset) |
| H | within-CBF HSC/stem-program EFS Kaplan–Meier | `panelG_cbf_stem_km.R` | `TARGET_bulk_program_scores.csv` |
| — | composite | `composite_fig2_v3.py` | → `Figure_2__leukemic_states.png` |

## Upstream compute (produces the prognosis + all-states table)
- `compute_wholecohort_prognosis.R` — per-state whole-cohort EFS Cox (BH-FDR<0.10 → 6 favorable, 8 poor)
  → `LS_wholecohort_cox_efs.csv`, `LS_wholecohort_prognosis.csv` (drives panels B–E).
- `km_table_all_states.R` — KM log-rank + continuous Cox for all 49 leukemic + 7 normal states,
  whole-cohort and CBF, EFS and OS → `KM_all_states_table.csv` (also feeds the Figure 3 volcano).
- CIBERSORTx (web, S-mode, QN off, relative, 100 perms; Jobs 24–27) is off-box; fractions in `LS_survival.csv`.

## Key corrected results (see [[project_v2_recluster_rebuild]])
- 14 states EFS-prognostic (6 favorable = CBF/APL/CEBPA GMP; 8 poor = proliferation/MLL). The old
  "LS12 stem-like drives relapse" story does NOT reproduce.
- Stemness does NOT separate prognosis at the state level (panel D, NS) — honest.
- Normal-HSC content is not prognostic (`NORM_HSPC`); the old "normal HSC prognostic" was a
  contamination artifact (old LS42). Within CBF, stem-program-high patients relapse more (panel H).

## Provenance
No fabricated values; every panel derives from the files above. Stale 48-state (contaminated) scripts
are archived in `legacy_v1_48state/`; the earlier PAC scheme in `legacy_PAC/`.
