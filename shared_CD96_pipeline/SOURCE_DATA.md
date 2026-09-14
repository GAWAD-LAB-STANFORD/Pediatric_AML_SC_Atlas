# Source data & reproducibility index

Per-figure mapping of **generating script → source-data file(s) → reproducible now?**
Audited against the package on disk. "✓ now" = renders/recomputes from files already in the
package; "⏳ handoff" = needs an item still being retrieved (see `MISSING_AND_EXTRAS.md`).

## Main figures
| Fig | Script | Source data (in package) | Status |
|----|--------|--------------------------|--------|
| **1** atlas | `scripts/legacy/Figure1/*`, `rebuild_figure1_umaps.py` | UMAPs ← `count_matrix/` (obs + X_umap); `source_data/legacy/Figure1/pca_…repca.csv`, MRD xlsx | UMAPs **✓ now**; dendrogram/patient-heatmap **⏳** need `AML-sample-infor.csv` (Table S1) + `AML_scRNA_Cell_summary_Cell_Type.csv` |
| **2** subclusters/PAC/GO | `scripts/legacy/Figure2/Figure2_clusterprofiler.R` | `source_data/legacy/Figure3/AML_PAC_Rank_gene_wilcoxon{,_logFC,_padj}.csv`, `…/Figure2/AML_Module_hub_gene_list.csv`, `AML_scRNA_Combo_PAC_cell_counts.csv`; UMAPs ← `count_matrix/` | **✓ now** (GO needs online msigdbr) |
| **3** TF regulons | `scripts/legacy/server/Scanpy_in_R.R` | regulon-AUC matrix `reg-py_3.8.csv` | **⏳ handoff** (on Oak) |
| **4** T/NK | `scripts/legacy/Figure4/*` | `…/Figure4/AML_scRNA_T_Cell_Cell_Type.csv`, `…_PAC_anno_by_Sample.csv`, Tcell rank tables; UMAPs ← `count_matrix/` | composition/markers **✓ now**; survival **⏳** (CIBERSORTx_T_NK results + TARGET clinical) |
| **5** target discovery (12 antigens) | `scripts/fig5.R` | `source_data/fig5*` (funnel, UMAP, head-to-head, cyto/PAC, TARGET + BeatAML bulk) | **✓ now** |
| **6** toxicity disqualification | `fig6_toxicity.R` + `fig6_DE.R` → `build_fig6_assemble.py` | `source_data/fig6tox_{organ,hema,deadly,celltypes}.csv`, `fig6D_radar.csv` (← `external_data/census/` + `count_matrix/`) | **✓ now** |
| **7** combination strategy | `cellchat/fig7_circos_3panel.R` + `fig7_coverage.R` + `fig7_combo.R` → `build_fig7_assemble.py` | `source_data/fig7_combo.csv`, `fig5E_coverage.csv`, `scripts/cellchat/cellchat_3path_input.csv` | **✓ now** |

## Supplementary figures
| Fig | Source data | Status |
|----|-------------|--------|
| **S1** healthy BM | `count_matrix/` | ✓ now |
| **S2** atlas covariates | `count_matrix/` | ✓ now |
| **S3** deconvolution perf. | DWLS/CIBERSORTx pseudobulk results | ⏳ (regenerable from counts + signature) |
| **S4** cluster survival HR | TARGET clinical + cluster fractions | ⏳ (TARGET clinical, public) |
| **S5** regulon detail | `reg-py_*` AUCell | ⏳ (Oak) |
| **S6** hdWGCNA (was Fig 4) | `source_data/legacy/FigS_hdWGCNA/*MEs*.tsv`, `PAC_anno_Module_data.rds` | **✓ now** |
| **S7** hdWGCNA detail | `…/FigS_hdWGCNA/*` | **✓ now** |
| **S8** composition / cycle | `count_matrix/` obs | ✓ now |
| **S9** T/NK detail + HR | CITE ADT + CIBERSORTx survival | ⏳ |
| ~~**S10** surfaceome MCC-F1~~ | **DROPPED (2026-08-27)** — replaced by the strict-LSC analysis (current Fig S8). Legacy code archived to `Archive_superseded/dropped_legacy_scripts_20260827/FigS_surfaceome/`; its MCCF1/CSPA source data to `Archive_superseded/dropped_MCC_screen_Figure_S8_source_data_20260827/`. | superseded |
| **S11** CITE-seq norm | CITE-seq ADT object (`…CITEseq…h5mu`) | ⏳ (Oak) |
| **S12** CD96 window violins | `source_data/fig6F_*` | **✓ now** |
| **S13** CD96 metacell corr | `source_data/figS7_*` | **✓ now** |
| **S14** CD96 CellChat | `source_data/cd96_*` | **✓ now** |
| **S15** Therapeutic window per PAC | `source_data/figS_pac_window.csv` (← `count_matrix/` via `build_pac_window.py`); `figS_pac_window.R` | **✓ now** |
| **S16** Therapeutic window per cytogenetic subtype | `source_data/figS_cyto_window.csv` (← `count_matrix/` via `build_cyto_window.py`); `figS_cyto_window.R` | **✓ now** |
| **S17** PPAC-specific discovery (% positive, tox≤10%) | `source_data/funnel_ppac_specific.csv` + `…_composition.csv` (← `count_matrix/` + `ref/surfaceome_genes.txt` via `build_funnel_ppac_specific.py`); `figS_ppac_specific.R` | **✓ now** |
| **S18** Per-PPAC efficacy-vs-toxicity (5C style) | `source_data/figS_ppac_headtohead.csv` (← `count_matrix/` + S17 ranking via `build_ppac_headtohead.py`); `figS_ppac_headtohead.R` | **✓ now** |
| **S19** Additive CD96 combinations per PPAC (6D style) | `source_data/figS_ppac_additive.csv` (← `count_matrix/` + S17 ranking via `build_ppac_additive_window.py`); `figS_ppac_additive.R` | **✓ now** |
| **Table S1** clinical | `AML-sample-infor.csv` | ⏳ (= the table itself) |

## Answers
1. **All figures reproducible from input + scripts?** The **CD96 figures (5–7, S12–S14) yes,
   fully** (0 missing inputs, no external paths). **Figs 1, 2, S1, S2, S6, S7, S8 also reproduce
   now** from the bundled `count_matrix/` + `source_data/legacy/`. The remainder (3, 4-survival,
   S3–S5, S9–S11) need the short Oak/Yakun handoff list and script re-pathing.
2. **Source-data files saved per figure?** Yes for every CD96 figure (one CSV per panel set in
   `source_data/`) and for most legacy figures (`source_data/legacy/<figure>/`). The few gaps
   are exactly the handoff items above.
3. **Datasets to add as Supplementary Data/Tables** (formalize these for submission):
   - **SD1 — raw count matrix** (`count_matrix/ALSF_AML_raw_counts.h5ad`) — primary; deposit in GEO.
   - **Table S1 — patient/sample clinical** (`AML-sample-infor.csv`).
   - **SD — cluster/PAC marker genes** (`AML_PAC_Rank_gene_wilcoxon*`).
   - **SD — hdWGCNA modules** (membership/hub + `*MEs*`).
   - **SD — SCENIC regulons** (regulon→target + AUCell; from Oak).
   - **SD — surfaceome MCC-F1 screen** (full ranked list).
   - **SD — CD96 & 12-antigen metrics** (`fig5*`, `fig6D_radar`, `fig7_combo`, `fig6E_coverage`).
   - **SD — CellChat CD96 communication** (`cd96_persample_long`, `…aggregated`; Fig S14).
   - **SD — body-wide normal expression (Census)** (`external_data/census/`, `fig6tox_*`).
4. **Base count matrix?** **Yes — it's the most important supplementary dataset and is already
   bundled** (`count_matrix/`, ~458 MB h5ad + `cell_metadata.csv.gz`); deposit the equivalent
   10x matrix in **GEO** with the cell metadata. It alone regenerates all scRNA panels.
