# Figure 3 — CBF-AML deep-dive (leukemic states refine event-free survival within CBF)

Data-parametrized pipeline. Everything reads paths from `config.R` / `config.py`, so
re-running on the corrected/new data is just pointing `FIG3_DATA` at the right directory.
Nothing is hardcoded to a specific state number — the **focus state** (top within-CBF
poor predictor) and the CBF-prognostic state list are chosen **from the data** each run.

## Panels
| panel | content | script |
|---|---|---|
| A | CBF-prognostic states × CBF cytogenetic enrichment + stemness + within-CBF EFS Cox HR(p) | `02_panelA_cbf_states.R` |
| B | Focus-state deconvolved-fraction tertiles → Kaplan–Meier of EFS | `03_panelB_focus_km.R` |
| C | Focus state as a distinct, targetable population (identity signature + surface antigens vs normal HSPC / favorable / poor) | `06_build_panelCD.py` → `07_panelsCD.R` |
| D | HOX gene expression across CBF-prognostic states | `06_build_panelCD.py` → `07_panelsCD.R` |
| E | GO:BP over-representation of per-state up-genes (CBF-prognostic states) | `04_compute_state_DE.py` → `05_panelE_cbf_go.R` |
| — | composite | `08_composite.py` |

## Run
```bash
bash run_all.sh                 # uses 08_Source_Data/Figure_3/v3
FIG3_DATA=/path/to/other bash run_all.sh   # swap datasets
```

## Inputs the data dir must contain
- `LS_survival.csv` — per-patient columns `LS_1..LS_K` (deconvolved fractions) + `sample`,
  `subtype2` (cytogenetic), and `os_*/efs_*` (or EFS is merged from the Figure-2 TARGET
  survival file). **Produced by stitching the CIBERSORTx output** (the one step done off-box).
- `LS_annotation.csv` — `LS, n_cells, lineage_program, stem_z, top_cyto, ...`
- `meancpm.npy`, `genes.npy`, `reference_assignments.csv`, `jaccard_percluster.csv`
- (derived, written by the pipeline) `LS_cbf_order.csv`, `LS_focus_state.txt`,
  `LS_upgenes_wilcoxon.csv`, `panelC_targetable.csv`, `panelD_hox.csv`

## Provenance / correctness notes
- Leukemic cells = Cell Type ∈ {AML, AML-MKI67, AML-PCNA, AML-CD1C} (the corrected set;
  the earlier bug included normal cells via `AML Sub-Clusters`).
- States = res-4.5 clusters with per-cluster Jaccard ≥ 0.35 (`config.py: stable_raws()`).
- No fabricated values; every panel derives from the files above.
