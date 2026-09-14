# Fig 4 (orig Fig 5) — T/NK-cell states + prognosis

Legacy figure carried from the original manuscript analysis (`Older_Files/ALSF AML scRNA-seq Paper/ALSF_AML_plot/`).

`Cell_Fraction_plot_summary.R` builds the T/NK fractions, marker plots and GZMK/GZMB-CD8 survival.

## TCR repertoire — Supplementary Figure S12

Built from the real Cell Ranger VDJ object (`ALSF_AML_plot/H5AD/ALSF_AML_total_Tcell_rna_withTCR.h5ad`):

1. `compute_tnk_tcr2.py` → writes the TCR source CSVs to `08_Source_Data/Figure_4/v2/`
   (`tcr_sharing.csv`, `tcr_exhaustion_by_expansion.csv`, `tcr_clonality_by_cyto.csv`,
   `tcr_clonality_downsampled.csv`, `tcr_trbv_usage.csv`).
2. `panel_tnk_tcr_extra.R` → renders the panels to `04_Main_Figures/Figure_4_panels/` (PDF + PNG).
3. `composite_supp_tcr.py` → composites the 3 deployed panels into
   `05_Supplementary_Figures/Figure_S12__TNK_TCR_repertoire.{pdf,png}` (run step 2 first).

Deployed S12 = **A** clonotype sharing (Jaccard) · **B** exhaustion vs clonal expansion (null) ·
**C** clonality by cytogenetics (null). The per-subset downsampled clonality (main Fig 4F analog) and
the TRBV-usage panel are built by step 2 but are **not** in the deployed supplement — TRBV was dropped
because the HBM normal reference is only 2 donors. Interpretation lives in the S12 legend, not on the panels.
