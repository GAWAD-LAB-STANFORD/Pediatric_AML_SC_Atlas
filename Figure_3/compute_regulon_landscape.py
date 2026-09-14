#!/usr/bin/env python
"""Regulon landscape data: re-aggregate the pySCENIC per-cell AUCell (integrated loom) by
the 49 corrected leukemic states AND the normal cell-type compartments (reference), then
z-score each regulon across all columns. SCENIC is label-independent, so re-aggregation by
new labels is equivalent to re-running. Writes the per-column z matrix. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
LOOM = os.path.join(os.path.dirname(C.H5), "Total_30_scenic_integrated-output_python_3.8_new.loom")
D3 = C.DATA
f = h5py.File(LOOM, "r")
cid = np.array([c.decode() if isinstance(c, bytes) else c for c in f["col_attrs"]["CellID"][:]])
ct  = np.array([c.decode() if isinstance(c, bytes) else c for c in f["col_attrs"]["Cell_Type"][:]])
auc = f["col_attrs"]["RegulonsAUC"][:]; regs = list(auc.dtype.names)
A = pd.DataFrame({r: auc[r] for r in regs}, index=cid); f.close()
# leukemic states
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}; ra["LS"] = ra["res4.50"].map(raw2ls)
bc2ls = dict(zip(ra.barcode, ra.LS))
# normal compartments from the loom Cell_Type labels (non-leukemic reference)
NORM = {"HSPC": "NORM_HSPC", "Myeloid Pro": "NORM_MyeloidPro",
        "CD14_Monocyte": "NORM_Monocyte", "CD16_Monocyte": "NORM_Monocyte", "Macrophage": "NORM_Macrophage",
        "mDC": "NORM_DC", "pDC": "NORM_DC", "Erythrocytes": "NORM_Erythroid",
        "Naive CD4T": "NORM_T_NK", "Naive CD8T": "NORM_T_NK", "CTL": "NORM_T_NK", "NK": "NORM_T_NK",
        "GZMB+CD8T": "NORM_T_NK", "MEM_CD8T": "NORM_T_NK",
        "CD20+B": "NORM_B", "ProB": "NORM_B", "PreB": "NORM_B", "PlasmaB": "NORM_B", "CD34+ProB": "NORM_B"}
col = []
for b, c in zip(A.index, ct):
    if b in bc2ls: col.append(bc2ls[b])
    elif c in NORM: col.append(NORM[c])
    else: col.append(None)
A["col"] = col; A = A[A["col"].notna()]
mean = A.groupby("col")[regs].mean()
z = (mean - mean.mean(0)) / (mean.std(0) + 1e-9)
z.to_csv(os.path.join(D3, "regulon_landscape_z.csv")); mean.to_csv(os.path.join(D3, "regulon_landscape_mean.csv"))
n_ls = sum(1 for c in z.index if c.startswith("LS_")); n_no = sum(1 for c in z.index if c.startswith("NORM_"))
print(f"wrote regulon_landscape_z.csv: {z.shape[1]} regulons x {z.shape[0]} columns ({n_ls} leukemic states + {n_no} normal compartments)")
