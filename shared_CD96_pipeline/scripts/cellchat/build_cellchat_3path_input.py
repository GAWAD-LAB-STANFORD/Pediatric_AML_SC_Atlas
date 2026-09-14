#!/usr/bin/env python
# ====================================================================
# build_cellchat_3path_input.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_cellchat_3path_input.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Build the CellChat input for Figure 7's three circos — the 3 selected targets that ARE
curated ligand/receptors in CellChatDB.human:
  CD96   : CD96 (lig, on AML) -> PVR / NECTIN1            (AML = sender)
  TNFRSF4: TNFSF4 (OX40L) -> TNFRSF4 (OX40, on AML)        (AML = receiver)
  ITGAX  : C3 / ICAM1 / FCER2A / THY1 -> ITGAX+ITGB2       (AML = receiver)
Per cell: compartment (same 8 hematopoietic groups as the CD96 CellChat), sample, raw counts.
-> scripts/cellchat/cellchat_3path_input.csv
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CM = os.path.join(ROOT, "count_matrix", "ALSF_AML_raw_counts.h5ad")
OUT = os.path.join(HERE, "cellchat_3path_input.csv")

WANT = ["CD96", "NECTIN1", "PVR", "NECTIN2", "PVRL2", "TNFSF4", "TNFRSF4", "C3",
        "ITGAX", "ITGB2", "FCER2A", "FCER2", "ICAM1", "THY1"]   # FCER2A may be FCER2; PVRL2 = NECTIN2 alias
COMPART = {
    **{c: "AML (leukemic)" for c in ["AML", "AML-MKI67", "AML-PCNA"]},
    "0_HSPC": "HSPC", "Myeloid_Pro": "Myeloid_Pro",
    **{c: "Mature myeloid" for c in ["CD14_Monocyte", "CD16_Monocyte", "mDC", "pDC", "Macrophage", "AML-CD14", "AML-CD1C"]},
    **{c: "B-lineage" for c in ["CD20+B", "ProB", "PreB", "CD34+ProB", "PlasmaB", "AML-B"]},
    **{c: "T cells" for c in ["Naïve_CD4T", "Naïve_CD8T", "CTL", "Activated_CD4T", "AML-CD4T", "AML-CTL", "AML-Naïve_CD8T"]},
    **{c: "NK" for c in ["NK", "AML-NK"]},
    **{c: "Erythroid" for c in ["Erythrocytes", "AML-Ery"]},
}

A = ad.read_h5ad(CM)
genes = [g for g in WANT if g in A.var_names]
absent = [g for g in WANT if g not in A.var_names]
print(f"present: {genes}\nabsent : {absent}")
ct = A.obs["Cell Type"].astype(str).values
comp = np.array([COMPART.get(c, "Erythroid") for c in ct])
samp = A.obs["Sample"].astype(str).values
X = A[:, genes].X
Xd = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)

df = pd.DataFrame({"cell": [f"c{i}" for i in range(A.n_obs)], "sample": samp, "compartment": comp})
for j, g in enumerate(genes):
    df[g] = Xd[:, j].astype(int)
df = df[~df["sample"].str.contains("Healthy", case=False)]
df.to_csv(OUT, index=False)
print(f"wrote {OUT}  ({len(df)} cells, {df['sample'].nunique()} AML samples, {len(genes)} genes)")
