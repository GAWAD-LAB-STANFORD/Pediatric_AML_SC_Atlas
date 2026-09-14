# Figure → script map (authoritative, deployed numbering)

Every **deployed** figure — main Figures 1–7 and Supplementary Figures **S1–S27** (the files in
`../04_Main_Figures/` and `../05_Supplementary_Figures/`) — maps to a **self-contained folder** here.
The supplement code folders and the `08_Source_Data/` folders now match the deployed S-numbers
(the earlier superseded/drifted numbering has been reconciled).

## Self-contained layout
Each CD96 figure folder and each supplement folder contains:
```
Figure_<N>/
  scripts/       renderer(s) + builder(s) + bundled shared deps (_style.R / config.py / manuscript_palette.R)
  source_data/   plot-level CSVs the renderer reads
  figures/       render output lands here (deployed copy lives in 04_/05_)
  README.md      what it is + how to reproduce
```
Renderers reproduce the figure from `source_data/` alone (`cd scripts && Rscript <renderer>.R`, writing
`../figures/`). Builders (`build_*.py`) regenerate `source_data/` from the raw single-cell atlas /
TARGET bulk (GEO/dbGaP; path in `config.py` `H5`) and therefore need the deposited raw data.

## Main figures
| Fig | Folder | Renderer(s) | Notes |
|-----|--------|-------------|-------|
| 1 | `Figure_1/` | `ALSF_AML_Dendrogram_plotting.R` (+ scanpy atlas UMAPs, upstream) | atlas + sample dendrogram |
| 2 | `Figure_2/` | `composite_fig2_v3.py` (+ `panel*`/`compute*`) | 49 leukemic states + prognosis |
| 3 | `Figure_3/` | `composite_fig3_ls10.py` (+ many `panel_*`/`compute_*`) | LS10 proliferation deep-dive |
| 4 | `Figure_4/` | `composite_fig4.py` (+ `panel_tnk*`/`annotate_compose_tnk.py`) | T/NK states + prognosis |
| 5 | `Figure_5/` | `fig5.R` (A–F) | CD96 target discovery |
| 6 | `Figure_6/` | `fig6_toxicity.R`, `fig6_DE.R` → `build_fig6_assemble.py` | marrow-sparing toxicity |
| 7 | `Figure_7/` | `cellchat/fig7_circos_3panel.R`, `fig7_combo.R`, `fig7d_bestcombo.R`, `fig7_allpairs.R` → `build_fig7_assemble.py` | combination strategy |

Figs 5–7 are also rebuildable end-to-end as one pipeline: `cd shared_CD96_pipeline && bash run_all.sh`.

## Supplementary figures (deployed S1–S27)
| Dep | Folder | Title | Renderer (in `scripts/`) | Provenance |
|-----|--------|-------|--------------------------|-----------|
| S1  | `Figure_S1/`  | Healthy bone-marrow reference | *external scanpy* | UMAP coords shipped; atlas to GEO |
| S2  | `Figure_S2/`  | AML atlas covariates | `panel_S2_composite.py` | reads atlas (GEO) |
| S3  | `Figure_S3/`  | Selection of the 49 reproducible leukemic states | `panel_state_selection.R` | in-tree |
| S4  | `Figure_S4/`  | Bulk-deconvolution validation | `panel_deconv_validation.R` | in-tree |
| S5  | `Figure_S5/`  | Transcription-factor regulon detail | *external pySCENIC* | Z-score matrices shipped; looms to GEO |
| S6  | `Figure_S6/`  | Cellular composition and cell cycle | `panel_S9_composite.R` | in-tree |
| S7  | `Figure_S7/`  | Annotation of the 49 leukemic states | `panel_state_annotation.R` | in-tree |
| S8  | `Figure_S8/`  | Whole-cohort EFS hazard ratios for the 49 states | `panel_state_hr_forest.R` | in-tree |
| S9  | `Figure_S9/`  | State-level stemness does not distinguish prognosis | `panelD_stem_by_prog.R` | in-tree |
| S10 | `Figure_S10/` | Absence of a primitive quiescent LSC population | `composite_figS_lsc.py` | panels read atlas (GEO) |
| S11 | `Figure_S11/` | CITE-seq normalization and gating | *external muon/ADT notebook* | plot data shipped; object to GEO |
| S12 | `Figure_S12/` | CITE-seq surface-protein validation of the T/NK subsets | `build_figS11_citeseq.py` | reads CITE-seq object (GEO) |
| S13 | `Figure_S13/` | T-cell receptor repertoire of the marrow T/NK subsets | `panel_tnk_tcr_extra.R` → `composite_supp_tcr.py` | in-tree (VDJ object to GEO for rebuild) |
| S14 | `Figure_S14/` | Lead-target cell-type specificity | `figS_leads_celltype.R` | in-tree |
| S15 | `Figure_S15/` | Poor-prognosis-state (LS) surface-target discovery | `figS_ls_specific.R` | in-tree |
| S16 | `Figure_S16/` | Per-poor-prognosis-state (LS) efficacy vs toxicity | `figS_ls_headtohead.R` | in-tree |
| S17 | `Figure_S17/` | CD96–comparator co-expression (metacells + bulk) | `figS.R` (correlations) | in-tree |
| S18 | `Figure_S18/` | AML immunotherapy clinical-trial landscape | `figS_aml_trials.R` | in-tree |
| S19 | `Figure_S19/` | Discovered antigens versus antigens in clinical trials | `figS_venn.R` | in-tree |
| S20 | `Figure_S20/` | Therapeutic-window violins, KMT2Ar/MLL-excluded | `figS.R` (nonMLL) | in-tree |
| S21 | `Figure_S21/` | Candidate antigens by cytogenetic subtype, pediatric TARGET | `figS_cyto_violin.R` | in-tree |
| S22 | `Figure_S22/` | Candidate antigens by cytogenetic subtype, adult BeatAML | `figS_cyto_violin.R` | in-tree |
| S23 | `Figure_S23/` | Adult AML single-cell target validation | `fig_adult_scrna.R` | in-tree |
| S24 | `Figure_S24/` | CD96 cell–cell communication (CellChat) + nectin axis | `figS8.R` | in-tree |
| S25 | `Figure_S25/` | Additive therapeutic window per poor-prognosis LS | `figS_ls_window.R` | in-tree |
| S26 | `Figure_S26/` | Additive therapeutic window by cytogenetic subtype | `figS_cyto_window.R` | in-tree |
| S27 | `Figure_S27/` | Every two-target combination landscape | `fig7_allpairs.R` | in-tree (also builds Fig 7D/E) |

**External figures** (no in-tree renderer; plot-level source data shipped in the folder, raw objects to
GEO): S1 (scanpy atlas UMAPs), S5 (pySCENIC AUCell), S11 (muon/ADT normalization notebook).

**Reconciled in this pass:** the Fig-4 T/NK trio is now ordered normalization (S11) → validation (S12)
→ TCR (S13); the former all-pairs S28 is now S27 (gap closed); the dropped surfaceome screen (old S12)
and the merged per-PPAC window (old S27) are archived under `Archive_superseded/`. The two hdWGCNA
figures (non-reproducing) are not in the deployed set; their scripts are in
`Archive_superseded/old_code_folder_stragglers/`.

## Reproducing
- **Any single figure:** `cd Figure_<N>/scripts && Rscript <renderer>.R` (writes `../figures/`), reading
  the bundled `../source_data/`. See each folder's `README.md`.
- **Whole CD96 story at once:** `cd shared_CD96_pipeline && bash run_all.sh` (canonical integrated
  pipeline; the per-figure CD96 folders are self-contained copies of its scripts + source data).
- See `README_code_locations.md` and `../MANIFEST.md` for the full inventory.
