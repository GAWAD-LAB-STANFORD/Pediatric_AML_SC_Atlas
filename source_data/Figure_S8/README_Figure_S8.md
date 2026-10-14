# Figure S8 — absence of a primitive quiescent leukemic stem-cell population

Rendered by the scripts in `../../07_Code/Figure_S8/scripts/`:
`compute_strict_lsc.py` applies the strict LSC gate to the atlas and writes
`strict_LSC_funnel.csv` and `strict_LSC_barcodes.csv`; `panel_strictLSC_markers.R`,
`panel_strict_lsc_funnel.R` and `panel_strict_lsc_umap.py` draw the panels;
`composite_figS_lsc.py` assembles them.

The input the renderers need, `strictLSC_marker_profile.csv`, ships in
`../../07_Code/Figure_S8/source_data/`. `strict_LSC_funnel.csv` and
`strict_LSC_barcodes.csv` are derived outputs that `compute_strict_lsc.py` regenerates
from the atlas; they are not deposits. `normal_HSC_gated_barcodes.csv` ships under
`../Figure_2/v3/`.

NOTE ON THE OTHER FILES IN THIS FOLDER. The `MCCF1_*`, `ALSF_AML_MCCF_*`,
`ALSF_RNA-seq_*`, `CSPA_human_surface_protein.csv` and `interaction_input.csv` files
are a legacy surfaceome MCC/F1 screen carried over from the original manuscript
analysis. No Figure S8 script reads them, and they do not correspond to Figure S8 as
numbered in the submitted manuscript. They are retained for provenance only.
