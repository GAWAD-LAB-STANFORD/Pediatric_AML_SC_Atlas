#!/usr/bin/env python
"""Pseudotime ordering data: diffusion pseudotime (dpt, rooted at a normal HSPC cell -- verified
iroot Cell_Type=0_HSPC) per cell, assigned to the 49 leukemic states OR the 8 normal hematopoietic
compartments (reference anchors). Palantir pseudotime/entropy are NOT used (mis-oriented / saturated
in this object). Writes per-cell (group, dpt) for the ordering figure. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"); import config as C
f = h5py.File(C.H5, "r"); dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
obs = dec(f["obs"]["_index"][:])
ct = pd.Series(dec(f["obs"]["Cell_Type"]["categories"][:])[f["obs"]["Cell_Type"]["codes"][:]], index=obs)
dpt = pd.Series(f["obs"]["dpt_pseudotime"][:], index=obs); f.close()
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
NORM = {"0_HSPC":"NORM_HSPC","Myeloid_Pro":"NORM_MyeloidPro","CD14_Monocyte":"NORM_Monocyte","CD16_Monocyte":"NORM_Monocyte",
        "Macrophage":"NORM_Macrophage","mDC":"NORM_DC","pDC":"NORM_DC","Erythrocytes":"NORM_Erythroid",
        "Naïve_CD4T":"NORM_T_NK","Naïve_CD8T":"NORM_T_NK","CTL":"NORM_T_NK","NK":"NORM_T_NK","Activated_CD4T":"NORM_T_NK",
        "CD20+B":"NORM_B","ProB":"NORM_B","PreB":"NORM_B","PlasmaB":"NORM_B","CD34+ProB":"NORM_B"}
grp = []
for b in obs:
    v = bc2ls.get(b)
    if isinstance(v, str) and v.startswith("LS_"): grp.append(v)
    else: grp.append(NORM.get(ct.get(b, ""), None))
d = pd.DataFrame({"group": grp, "dpt": dpt.values}, index=obs).dropna(subset=["group"])
d.to_csv(os.path.join(C.DATA, "pseudotime_by_group.csv"))
print("cells:", len(d), "| leukemic:", d.group.str.startswith("LS_").sum(), "| normal:", d.group.str.startswith("NORM_").sum())
print("normal compartment median dpt (HSPC should be earliest):")
print(d[d.group.str.startswith("NORM_")].groupby("group")["dpt"].median().round(3).sort_values().to_string())
