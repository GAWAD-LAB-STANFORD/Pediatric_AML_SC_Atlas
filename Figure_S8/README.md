# Figure S8 — Absence of a primitive quiescent LSC population

**Self-contained.** Reproduces deployed **Figure S10** (`../../05_Supplementary_Figures/Figure_S8__LSC_scarcity.pdf`) from its own source data.
Part of the atlas/prognosis analysis stream (legacy Figure_1-4 pipelines); scripts here are copies so this figure stands alone.

## Reproduce
```bash
cd scripts
python composite_figS_lsc.py
```
Reads `../source_data/` and writes `../figures/`.

## Note
Panel builders read the scanpy atlas (config H5); composite stitches the rendered panels.

## Contents
- `scripts/` — composite_figS_lsc.py, compute_strict_lsc.py, panel_strictLSC_markers.R, panel_strict_lsc_funnel.R, panel_strict_lsc_umap.py + shared manuscript_palette.R/config.py
- `source_data/` — 15 file(s): ALSF_AML_MCCF_AML_AML_0.2.csv, ALSF_AML_MCCF_Prognosis_vs_Gene_AML_0.2.csv, ALSF_RNA-seq_HSPC_low_PACs_0.1.csv, ALSF_RNA-seq_HSPC_low_other_0.1.csv, ALSF_RNA-seq_mccf1_specific_gene_NK.csv, ALSF_RNA-seq_mccf1_specific_gene_PlasmaB.csv, ALSF_RNA-seq_mccf1_specific_gene_pDC.csv, CSPA_human_surface_protein.csv, MCCF1_AML_specific_gene_list copy.csv, MCCF1_AML_specific_gene_list.csv, MCCF1_AML_specific_gene_list_NCI_PAML_validation.csv, MCCF1_mcc_0.5_Cell_surface_marker_full_list.csv, MCCF1_mcc_0.5_Cell_surface_marker_full_list.txt, README_Figure_S8.md, interaction_input.csv
