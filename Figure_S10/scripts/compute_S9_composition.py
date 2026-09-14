#!/usr/bin/env python
"""S9 (composition + cell cycle) DATA, rebuilt from the current h5ad obs. Per-sample lineage
composition, per-sample cell-cycle phase fractions with outcome, and normal-lineage abundance
(AML vs healthy BM). Cell-cycle-by-prognosis (Panel D) is computed separately
(state_cellcycle_by_prognosis.csv). No fabricated values."""
import os, numpy as np, pandas as pd, scanpy as sc
H5 = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
D2 = "../source_data"
a = sc.read_h5ad(H5, backed="r")
o = a.obs[["SampleID","lineage","phase","Prognosis","SampleType"]].astype(str).copy()

# A: per-sample lineage composition
lin = (o.groupby("SampleID")["lineage"].value_counts(normalize=True).rename("frac").reset_index())
lin.to_csv(os.path.join(D2,"s9_lineage_by_sample.csv"), index=False)

# B/C: per-sample cell-cycle phase fractions + outcome
ph = o.groupby("SampleID")["phase"].value_counts(normalize=True).unstack(fill_value=0)
for c in ["G1","S","G2M"]:
    if c not in ph: ph[c]=0.0
outc = o.groupby("SampleID")["Prognosis"].first()
st   = o.groupby("SampleID")["SampleType"].first()
ph["outcome"] = outc.reindex(ph.index).values
ph["sample_type"] = st.reindex(ph.index).values
ph = ph.reset_index()
ph["outcome"] = ph["outcome"].replace({"0_HealthyBM":"HBM"})
ph.to_csv(os.path.join(D2,"s9_phase_by_sample.csv"), index=False)

# E: normal-lineage abundance AML vs HBM (depletion)
NORMAL = ["0_HSPC","Myeloid_Pro","Erythrocytes","Monocyte","NK","T","B","PlasmaB"]
lin2 = lin.copy(); lin2["sample_type"] = st.reindex(lin2["SampleID"]).values
lin2["sample_type"] = lin2["sample_type"].replace({"HealthyBM":"HBM"})
norm = lin2[lin2["lineage"].isin(NORMAL)].copy()
# ensure zero rows for samples lacking a normal lineage
full = pd.MultiIndex.from_product([lin2["SampleID"].unique(), NORMAL], names=["SampleID","lineage"]).to_frame(index=False)
norm = full.merge(norm, on=["SampleID","lineage"], how="left")
norm["frac"] = norm["frac"].fillna(0.0)
norm["sample_type"] = st.reindex(norm["SampleID"]).values
norm["sample_type"] = norm["sample_type"].replace({"HealthyBM":"HBM"})
norm.to_csv(os.path.join(D2,"s9_normal_lineage_by_sample.csv"), index=False)
print("samples:", ph["SampleID"].nunique(), "| AML:", (st=='AML').sum(), "HBM:", (st=='HealthyBM').sum())
print("outcome per sample:", ph["outcome"].value_counts().to_dict())
print("wrote s9_lineage_by_sample, s9_phase_by_sample, s9_normal_lineage_by_sample")
