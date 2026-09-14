#!/usr/bin/env python
# ====================================================================
# build_fig6A_spec_leads.py  |  Figure 6 (builder)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_fig6A_spec_leads.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""ALTERNATE Figure 6 panel A (cell-type specificity) for the 5 LEAD targets.
% positive + mean expression (CP10k log1p) of CD96/CD9/SUCNR1/AMN/IL1RAP across the
8 hematopoietic compartments. -> source_data/fig7B_spec.csv (overwrites; 5 leads).
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
CM = os.path.join(ROOT, "count_matrix", "ALSF_AML_raw_counts.h5ad")
OUT = os.path.join(ROOT, "source_data", "fig7B_spec.csv")
LEADS = ["CD96", "CD9", "SUCNR1", "AMN", "IL1RAP"]

COMPART = {  # compartment -> cell types (canonical only; matches the published panel)
    "AML (leukemic)": ["AML", "AML-MKI67", "AML-PCNA"],
    "HSPC": ["0_HSPC"],
    "Myeloid_Pro": ["Myeloid_Pro"],
    "Mature myeloid": ["CD14_Monocyte", "CD16_Monocyte", "mDC", "pDC", "Macrophage"],
    "B-lineage": ["CD34+ProB", "ProB", "PreB", "CD20+B", "PlasmaB"],
    "T cells": ["Naïve_CD4T", "Naïve_CD8T", "CTL", "Activated_CD4T"],
    "NK": ["NK"],
    "Erythroid": ["Erythrocytes"],
}

A = ad.read_h5ad(CM)
ct = A.obs["Cell Type"].astype(str).values
tot = np.asarray(A.X.sum(1)).ravel().astype(float); tot[tot == 0] = 1.0
Xl = A[:, LEADS].X
Ld = np.asarray(Xl.todense()) if sp.issparse(Xl) else np.asarray(Xl)
Lpos = Ld > 0
Lnorm = np.log1p(Ld / tot[:, None] * 1e4)

rows = []
for comp, cells in COMPART.items():
    m = np.isin(ct, cells)
    if m.sum() < 20:
        continue
    for i, g in enumerate(LEADS):
        rows.append(dict(gene=g, compartment=comp, pct=round(100 * Lpos[m][:, i].mean(), 3),
                         label=g, mean_expr=round(float(Lnorm[m][:, i].mean()), 5)))
df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(f"wrote {OUT}")
print(df.pivot(index="compartment", columns="gene", values="pct").reindex(columns=LEADS).to_string())
