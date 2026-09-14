#!/usr/bin/env python
# ====================================================================
# build_fig5D_panel.py  |  Figure 5 (builder)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_fig5D_panel.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""ALTERNATE Figure 5 panel D — coverage of the top-10 composite candidates + 6 clinical AML
antigens + FLT3 (control) across cytogenetic subtypes and PAC clusters.
% positive + mean expression (CP10k log1p): cyto over core-AML cells, PAC over all cells.
-> source_data/fig5D_cyto.csv, fig5D_pac.csv
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
CM = os.path.join(ROOT, "count_matrix", "ALSF_AML_raw_counts.h5ad")
SD = os.path.join(ROOT, "source_data")

# symbol -> display label.  composite top-10, then 6 clinical, then FLT3 control,
# then 3 PPAC-specific markers from S17 (TNFRSF4/ABCA7/ITGAX; CD7 already above).
GENES = {"CD96": "CD96", "CD9": "CD9", "SUCNR1": "SUCNR1", "AMN": "AMN", "IL1RAP": "IL1RAP",
         "LMBR1L": "LMBR1L", "NRXN2": "NRXN2", "FURIN": "FURIN", "UMODL1": "UMODL1", "CD7": "CD7",
         "CD33": "CD33", "IL3RA": "CD123", "CLEC12A": "CLL-1", "CD38": "CD38", "CD70": "CD70",
         "FLT3": "FLT3", "TNFRSF4": "TNFRSF4", "ABCA7": "ABCA7", "ITGAX": "ITGAX"}
CORE = ["AML", "AML-MKI67", "AML-PCNA"]
PACORD = ["0_HSPC", "Myeloid_Pro", "FPAC_1", "FPAC_2", "PPAC_1", "PPAC_2", "PPAC_3", "PPAC_4", "PPAC_5"]
CYTO_DROP = {"0_HealthyBM1", "0_HealthyBM2"}

A = ad.read_h5ad(CM)
ct = A.obs["Cell Type"].astype(str).values
cyto = A.obs["Cytogenetic"].astype(str).values
pac = A.obs["Prognosis-Associated Clusters"].astype(str).values
core = np.isin(ct, CORE)
tot = np.asarray(A.X.sum(1)).ravel().astype(float); tot[tot == 0] = 1.0
syms = list(GENES)
X = A[:, syms].X
Xd = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
Lpos = Xd > 0
Lnorm = np.log1p(Xd / tot[:, None] * 1e4)

def coverage(groups, mask_of, subset=None):
    rows = []
    for g in groups:
        m = mask_of(g)
        if subset is not None:
            m = m & subset
        if m.sum() < 20:
            continue
        for i, sym in enumerate(syms):
            rows.append(dict(gene=sym, group=g, pct=round(100 * Lpos[m][:, i].mean(), 2),
                             label=GENES[sym], mean_expr=round(float(Lnorm[m][:, i].mean()), 4)))
    return pd.DataFrame(rows)

cyto_groups = [c for c in pd.unique(cyto) if c not in CYTO_DROP]
coverage(cyto_groups, lambda g: cyto == g, subset=core).to_csv(os.path.join(SD, "fig5D_cyto.csv"), index=False)
coverage(PACORD, lambda g: pac == g).to_csv(os.path.join(SD, "fig5D_pac.csv"), index=False)
print(f"wrote fig5D_cyto.csv ({len(cyto_groups)} subtypes) + fig5D_pac.csv ({len(PACORD)} clusters), "
      f"{len(GENES)} genes (10 composite + 6 clinical + FLT3)")
