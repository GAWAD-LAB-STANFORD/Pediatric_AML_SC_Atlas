Superseded PPAC-scheme analysis. NOT part of the submitted figure set.

These files implement the prognosis-associated-cluster (PPAC/PAC) scheme that the 49
leukemic states (LS) replaced. They are retained for provenance only and are excluded
from 13_Source_Data.zip and 14_Code.zip.

Current equivalents, which do ship:
  funnel_ppac_specific.csv      -> funnel_ls_specific.csv     (Figure S16)
  figS_ppac_headtohead.csv      -> figS_ls_headtohead.csv     (Figure S15)
  figS_pac_window.csv           -> figS_ls_window.csv         (Figure S25)
  figS_ppac_additive.csv        -> no current equivalent; the analysis was retired
  build_funnel_ppac_specific.py -> build_funnel_ls_specific.py

SOURCE_DATA.md is also here. It mapped source data to an older supplementary numbering
(its "S15/S17/S18/S19" are not the submitted S15-S19) and referenced count_matrix/,
source_data/legacy/ and MISSING_AND_EXTRAS.md, none of which exist in this package.
FIGURE_SCRIPT_MAP.md is the authoritative map for the deployed numbering.

Also note: fig5D_pac.csv under this scheme had 81 rows; the shipped
08_Source_Data/Figure_5/fig5D_pac.csv has 171 rows and is the current file.
