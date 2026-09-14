#!/usr/bin/env python
# ====================================================================
# build_cellchat_IL1_input.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_cellchat_IL1_input.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Build the IL-1-pathway CellChat input for the alternate Figure 6 (IL1RAP panel).
Per cell: compartment (8 hematopoietic groups, same mapping as the CD96 CellChat),
sample, and raw counts of the IL-1 pathway genes present in the atlas. AML is the
RECEIVER here (IL1RAP is a receptor subunit). -> scripts/cellchat/cellchat_IL1_input.csv
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CM = os.path.join(ROOT, "count_matrix", "ALSF_AML_raw_counts.h5ad")
OUT = os.path.join(HERE, "cellchat_IL1_input.csv")

# IL-1 pathway genes present in the atlas (ligands + receptor subunits)
GENES = ["IL1A", "IL1B", "IL36A", "IL36G", "IL1R1", "IL1RL1", "IL1RL2", "IL1RAP"]
COMPART = {
    **{c: "AML (leukemic)" for c in ["AML", "AML-MKI67", "AML-PCNA"]},
    **{c: "HSPC" for c in ["0_HSPC"]},
    **{c: "Myeloid_Pro" for c in ["Myeloid_Pro"]},
    **{c: "Mature myeloid" for c in ["CD14_Monocyte", "CD16_Monocyte", "mDC", "pDC", "Macrophage", "AML-CD14", "AML-CD1C"]},
    **{c: "B-lineage" for c in ["CD20+B", "ProB", "PreB", "CD34+ProB", "PlasmaB", "AML-B"]},
    **{c: "T cells" for c in ["Naïve_CD4T", "Naïve_CD8T", "CTL", "Activated_CD4T", "AML-CD4T", "AML-CTL", "AML-Naïve_CD8T"]},
    **{c: "NK" for c in ["NK", "AML-NK"]},
    **{c: "Erythroid" for c in ["Erythrocytes", "AML-Ery"]},
}

A = ad.read_h5ad(CM)
ct = A.obs["Cell Type"].astype(str).values
comp = np.array([COMPART.get(c, "Erythroid") for c in ct])
samp = A.obs["Sample"].astype(str).values
X = A[:, GENES].X
Xd = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)

df = pd.DataFrame({"cell": [f"c{i}" for i in range(A.n_obs)], "sample": samp, "compartment": comp})
for j, g in enumerate(GENES):
    df[g] = Xd[:, j].astype(int)
# drop the 2 healthy-BM samples (CellChat is within AML samples, matching the CD96 run)
df = df[~df["sample"].str.contains("Healthy", case=False)]
df.to_csv(OUT, index=False)
print(f"wrote {OUT}  ({len(df)} cells, {df['sample'].nunique()} AML samples)")
print("  IL1RAP+ leukemic fraction per sample (targetable >=10%):")
aml = df[df.compartment == "AML (leukemic)"]
frac = aml.groupby("sample").apply(lambda x: 100 * (x["IL1RAP"] > 0).mean()).sort_values(ascending=False)
print(f"    targetable (>=10% IL1RAP+): {int((frac >= 10).sum())} of {len(frac)} samples")
print(frac.round(1).head(10).to_string())
