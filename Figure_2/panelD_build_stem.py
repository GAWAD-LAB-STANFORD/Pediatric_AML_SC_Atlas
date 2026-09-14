#!/usr/bin/env python
"""Figure 2D build (v3): extract the precomputed per-cell LSC (Eppert) and HSC scores
from the h5ad obs (LSC_EPPERT_Score, HSC_Score) for the corrected 49-state leukemic
cells, label each by its res-4.5 stable state, and write per-cell scores. No scores are
recomputed or fabricated -- they are read straight from the object."""
import os, sys, numpy as np, pandas as pd, h5py
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C
D2 = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"

# per-cell scores keyed by barcode
f = h5py.File(C.H5, "r")
bc = [x.decode() if isinstance(x, bytes) else x for x in f["obs"]["_index"][:]]
eppert = f["obs"]["LSC_EPPERT_Score"][:]; hsc = f["obs"]["HSC_Score"][:]; f.close()
sc = pd.DataFrame({"barcode": bc, "LSC_Eppert": eppert, "HSC": hsc}).set_index("barcode")

# corrected leukemic cells -> stable state
ra = pd.read_csv(C.REFA)
ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, lslabs = C.stable_raws(); raw2ls = {rw: lslabs[i] for i, rw in enumerate(raws)}
ra["LS"] = ra["res4.50"].map(raw2ls)
ra = ra[ra["LS"].notna()].copy()                       # cells in a stable state
ra = ra.join(sc, on="barcode")
miss = int(ra[["LSC_Eppert", "HSC"]].isna().any(axis=1).sum())
ra = ra.dropna(subset=["LSC_Eppert", "HSC"])
ra[["LS", "LSC_Eppert", "HSC"]].to_csv(os.path.join(D2, "LS_stem_percell.csv"), index=False)
med = ra.groupby("LS")[["LSC_Eppert", "HSC"]].median().reset_index()
med.to_csv(os.path.join(D2, "LS_stem_by_state.csv"), index=False)
print(f"wrote per-cell ({len(ra)} cells, {miss} dropped for missing score) and per-state ({len(med)} states) stem scores")
print(med.describe().round(3).to_string())
