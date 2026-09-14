#!/usr/bin/env python
"""Export raw-count matrix + metadata for a Monocle3 trajectory of the 49 leukemic states + 8 normal
compartments. Subsampled to <=1000 cells/group (small states kept whole) to keep it tractable while
preserving structure. Writes MTX (genes x cells), cells.csv (group, Cell_Type, SampleID batch, dpt,
existing UMAP) and genes.csv. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"); import config as C
OUT = "/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad/monocle_input"
os.makedirs(OUT, exist_ok=True)
f = h5py.File(C.H5, "r"); dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
obs = dec(f["obs"]["_index"][:]); var = dec(f["var"]["_index"][:])
ct = pd.Series(dec(f["obs"]["Cell_Type"]["categories"][:])[f["obs"]["Cell_Type"]["codes"][:]], index=obs)
sid = pd.Series(dec(f["obs"]["SampleID"]["categories"][:])[f["obs"]["SampleID"]["codes"][:]], index=obs)
dpt = pd.Series(f["obs"]["dpt_pseudotime"][:], index=obs)
um = f["obsm"]["X_umap"][:]
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
NORM = {"0_HSPC":"NORM_HSPC","Myeloid_Pro":"NORM_MyeloidPro","CD14_Monocyte":"NORM_Monocyte","CD16_Monocyte":"NORM_Monocyte",
        "Macrophage":"NORM_Macrophage","mDC":"NORM_DC","pDC":"NORM_DC","Erythrocytes":"NORM_Erythroid",
        "Naïve_CD4T":"NORM_T_NK","Naïve_CD8T":"NORM_T_NK","CTL":"NORM_T_NK","NK":"NORM_T_NK","Activated_CD4T":"NORM_T_NK",
        "CD20+B":"NORM_B","ProB":"NORM_B","PreB":"NORM_B","PlasmaB":"NORM_B","CD34+ProB":"NORM_B"}
grp = np.array([bc2ls.get(b) if (isinstance(bc2ls.get(b), str) and str(bc2ls.get(b)).startswith("LS_")) else NORM.get(ct.get(b,""), None) for b in obs])
meta = pd.DataFrame({"barcode": obs, "group": grp, "Cell_Type": ct.values, "SampleID": sid.values,
                     "dpt": dpt.values, "umap1": um[:,0], "umap2": um[:,1]})
meta = meta[meta.group.notna()].reset_index(drop=True)
rng = np.random.RandomState(0)
keep = meta.groupby("group").apply(lambda d: d.sample(min(len(d), 1000), random_state=0)).reset_index(drop=True)
keepset = set(keep.barcode)
mask = np.array([b in keepset for b in obs])
# raw counts
g = f["layers"]["Raw_Counts"] if "Raw_Counts" in f["layers"] else f["layers"]["counts"]
X = csr_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g.attrs["shape"]))
f.close()
Xs = X[mask]                                        # cells x genes
sel_obs = obs[mask]
gene_ok = np.asarray((Xs > 0).sum(0)).ravel() >= 10
Xs = Xs[:, gene_ok]; genes = var[gene_ok]
# reorder metadata to match sel_obs
keep = keep.set_index("barcode").loc[sel_obs].reset_index()
mmwrite(os.path.join(OUT, "counts.mtx"), Xs.T.tocoo())   # genes x cells
keep.to_csv(os.path.join(OUT, "cells.csv"), index=False)
pd.DataFrame({"gene": genes}).to_csv(os.path.join(OUT, "genes.csv"), index=False)
print("exported cells:", Xs.shape[0], "genes:", Xs.shape[1])
print("per-group counts:\n", keep.group.value_counts().sort_index().to_string())
