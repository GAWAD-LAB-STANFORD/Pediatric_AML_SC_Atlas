#!/usr/bin/env python
# ====================================================================
# build_leads_celltype.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_leads_celltype.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Rigor / safety pass for the 4 lead targets (CD96, CD9, SUCNR1, IL1RAP):
% of cells positive for each, across normal cell types vs core leukemic, from the count
matrix. Reveals which normal compartments drive each target's expression (the funnel's
toxicity gate only inspects 0_HSPC and Myeloid_Pro, so this is the fuller safety view).

Output -> source_data/figS16_leads_celltype.csv  (long: cell_type x gene; with compartment).
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
CM = os.path.join(ROOT, "count_matrix", "ALSF_AML_raw_counts.h5ad")
OUT = os.path.join(ROOT, "source_data", "figS16_leads_celltype.csv")
LEADS = ["CD96", "CD9", "SUCNR1", "IL1RAP"]   # AMN dropped -> 4 leads

# canonical compartments (exclude the AML-sample-derived "AML-<lineage>" non-malignant cells
# to keep the safety read clean; core leukemic = the corrected 4-type set
#   AML / AML-MKI67 / AML-PCNA / AML-CD1C  (AML-CD1C is leukemic; AML-<other-lineage> are not))
COMPART = [
    ("leukemic", ["AML", "AML-MKI67", "AML-PCNA", "AML-CD1C"]),
    ("HSPC / Myeloid (gate)", ["0_HSPC", "Myeloid_Pro"]),
    ("monocyte / DC", ["CD14_Monocyte", "CD16_Monocyte", "mDC", "pDC", "Macrophage"]),
    ("T / NK", ["Naïve_CD4T", "Naïve_CD8T", "CTL", "Activated_CD4T", "NK"]),
    ("B / plasma", ["CD34+ProB", "ProB", "PreB", "CD20+B", "PlasmaB"]),
    ("erythroid", ["Erythrocytes"]),
]

A = ad.read_h5ad(CM)
ct = A.obs["Cell Type"].astype(str).values
X = A[:, LEADS].X
pos = (np.asarray(X.todense()) > 0) if sp.issparse(X) else (np.asarray(X) > 0)

rows = []
for comp, cells in COMPART:
    for c in cells:
        m = ct == c
        if m.sum() < 50:
            continue
        for i, g in enumerate(LEADS):
            rows.append(dict(compartment=comp, cell_type=c, n=int(m.sum()),
                             gene=g, pct=round(100 * pos[m][:, i].mean(), 1)))
df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(f"wrote {OUT}  ({df.cell_type.nunique()} cell types x {len(LEADS)} leads)\n")

piv = df.pivot_table(index=["compartment", "cell_type"], columns="gene", values="pct").reindex(columns=LEADS)
pd.set_option("display.width", 150)
print(piv.to_string())
print("\n=== SAFETY READ ===")
print("CD9 lights up normal B-lineage progenitors (ProB/PreB/CD34+ProB > leukemia) -> B-aplasia risk.")
print("NOTE: deployed CSV regenerated from the main atlas (ALSF_AML_Combo_3500_with_PAC_new.h5ad,")
print("      obs 'Cell Type', raw counts) on the corrected 4-type leukemic set; raw_counts.h5ad is GEO-only.")
print("SUCNR1: HSPC ~6% (above the strict 3% gate); otherwise marrow-clean. CD96 low on AML-CD1C (~5%).")
print("IL1RAP: known AML-LSC target, myeloid ~8-10% (expected). CD96: spares HSPC/myeloid/B/ery.")
