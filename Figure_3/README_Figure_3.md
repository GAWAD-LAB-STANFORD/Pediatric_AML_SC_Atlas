# Figure 3 — LS10, an adverse proliferating myeloid-committed (GMP-like) leukemic state

Built on the corrected 49 leukemic states (Cell Type ∈ {AML, AML-MKI67, AML-PCNA, AML-CD1C};
Leiden res 4.5, per-cluster Jaccard ≥ 0.35). Shared config: `config.py`/`config.R`
(`stable_raws()` maps res-4.5 stable clusters → LS_1..LS_49 by size). Palette: `manuscript_palette.R`.
Data dir: `../../08_Source_Data/Figure_3/v3/`. Panels render to `../../04_Main_Figures/Figure_2_panels/`.

## Story
LS10 is the top whole-cohort EFS-poor state. It is **not** undifferentiated and **not** a quiescent
LSC: it is a **proliferating, myeloid-committed (MPO+, GMP-signature) blast that aberrantly retains
stem markers** (CD34/HOXA9 on its MPO+ cells at ~5× the normal-GMP rate). Its prognosis is confounded
with adverse cytogenetics. A strict, literature-grounded LSC search finds no primitive quiescent-LSC
population (supplementary figure).

## Main figure — panel → render script → compute → data
| panel | content | render | compute → data |
|---|---|---|---|
| A | selection volcano (all 49 states, whole-cohort EFS) | `panel_state_volcano.R` | `../Figure_2/km_table_all_states.R` → `Figure_2/v3/KM_all_states_table.csv` |
| B | LS10 vs NS-prolif (LS39/LS3) vs normal HSC, UMAP | `panel_ls10_umap.py` | reference_assignments + `Figure_2/v3/normal_HSC_gated_barcodes.csv` |
| C | differential GO:BP (LS10 vs LS39/LS3) | `panel_ls10_pathways.R` | `compute_ls10_de.py` → `LS10_up_nonhistone.csv`, `LS10_dn.csv` |
| D | LS10 whole-cohort EFS Kaplan–Meier | `km_state.R LS_10` | `LS_survival.csv` |
| E | LS10 deconvolved % across cytogenetic subtypes | `panel_ls10_by_cyto.R` | `LS_survival.csv` |
| F | HSC-sparing surface markers | `panel_ls10_surface.R` | `compute_ls10_surface.py` → `LS10_surface_dotplot.csv` |
| G | 4-module program map (quiescence/leuk-stem/prolif/diff) | `panel_ls_quiescence.R` | `compute_program_modules.py` → `LS_program_modules.csv` |
| — | composite | `composite_fig3_ls10.py` | → `Figure_3__LS10_proliferation.png` |

## Supplementary figure (LSC scarcity) — `../../05_Supplementary_Figures/Figure_S_LSC_scarcity.png`
| panel | content | render | compute → data |
|---|---|---|---|
| A | strict-LSC gating funnel (AML vs healthy BM) | `panel_strict_lsc_funnel.R` | `compute_strict_lsc.py` → `strict_LSC_funnel.csv` (counts also hardcoded in the R for the labels) |
| B | strict-LSC candidates on the UMAP | `panel_strict_lsc_umap.py` | `compute_strict_lsc.py` → `strict_LSC_barcodes.csv` |
| C | candidate marker profile vs HSC / GMP (hybrid, GMP-committed) | `panel_strictLSC_markers.R` | `compute_strict_lsc.py` → `strictLSC_marker_profile.csv` |
| — | composite | `composite_figS_lsc.py` | |

Strict true-LSC = CD34+CD38- AND quiescent (HLF/AVP+) AND non-cycling AND LSC-aberrant marker
(TIM3/CD96/CLEC12A, which normal HSC lack; Ding 2017; Kikushige PNAS 2010; CLEC12A/CD96 lit).

## Provenance / correctness
- Every panel derives from the files above; no fabricated values.
- LS10 poor-prognosis is CONFOUNDED with adverse cytogenetics (panel E), not significant within CBF,
  and rides on ~9% of TARGET patients — stated, not hidden.
- `legacy_cbf_deepdive/` = the earlier CBF Figure 3 (CBF now lives in Figure 2 panels F–H); kept for provenance.
