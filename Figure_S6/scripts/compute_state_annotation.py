#!/usr/bin/env python
"""S10 (leukemic-state annotation) DATA for the corrected 49-state scheme. Per-state mean
log-normalized expression of curated lineage-program markers (h5ad 'counts' layer), z-scored
per gene across the 49 states; plus per-state cross-cytogenetic breadth (effective # subtypes =
exp Shannon entropy of the per-cell Cytogenetic distribution) and a patient-dominance class.
Columns/states are defined exactly as the rest of Fig 2/3 (REFA res4.50 -> stable LS). Replaces
the stale legacy 48-state figure. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"); import config as C
D2 = C.DATA.replace("Figure_3", "Figure_2")

PROG = {  # curated hematopoietic lineage programs (Cell-cycle -> proliferative)
 "HSC/MPP":["CD34","HLF","AVP","CRHBP","HOXA9","MEIS1","MLLT3","PROM1","MECOM","GATA2"],
 "GMP":["MPO","ELANE","AZU1","PRTN3","CTSG","CEBPE","LYST","CFD"],
 "Monocyte":["LYZ","CD14","CSF1R","FCN1","VCAN","S100A8","S100A9","CD68"],
 "Macrophage":["C1QA","C1QB","C1QC","MRC1","SIGLEC1","TREM1"],
 "cDC":["CD1C","FCER1A","CLEC10A","CD207","CD1E","IRF8"],
 "pDC":["IRF7","LILRA4","GZMB","IL3RA"],
 "Megakaryocyte":["PPBP","PF4","ITGA2B","GP9","GP1BA","VWF","TUBB1","PF4V1"],
 "Erythroid":["GATA1","KLF1","HBB","HBA1","GYPA","ALAS2","TFRC","AHSP"],
 "Neutrophil":["FUT4","CEACAM8","MMP8","CD177","LTF","LCN2","CAMP","MMP9"],
 "B/plasma":["CD19","MS4A1","CD79A","CD79B","IGHM","JCHAIN","BLNK"],
 "T/NK":["CD3D","CD3E","IL7R","NKG7","GNLY","CD2","GZMH"],
 "proliferative":["MKI67","TOP2A","CDK1","UBE2C","CENPA","BIRC5","PCNA","NUSAP1","BUB1B"],
}
genes = [g for v in PROG.values() for g in v]

f = h5py.File(C.H5, "r")
def dec(a): return np.array([x.decode() if isinstance(x, bytes) else x for x in a])
var = dec(f["var"]["_index"][:]); obs = dec(f["obs"]["_index"][:])
def catcol(name):
    g = f["obs"][name]; return pd.Series(dec(g["categories"][:])[g["codes"][:]], index=obs)
cyto = catcol("Cytogenetic"); samp = catcol("SampleID")
gpos = {g: i for i, g in enumerate(var)}; present = [g for g in genes if g in gpos]
gi = [gpos[g] for g in present]
grp = f["layers"]["counts"]; shape = tuple(grp.attrs["shape"])
X = csr_matrix((grp["data"][:], grp["indices"][:], grp["indptr"][:]), shape=shape)
f.close()
E = pd.DataFrame(np.asarray(X[:, gi].todense()), index=obs, columns=present)

ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
E["LS"] = [bc2ls.get(b) for b in E.index]
lscol = pd.Series([bc2ls.get(b) for b in obs], index=obs)
E = E[E["LS"].apply(lambda v: isinstance(v, str) and v.startswith("LS_"))]

# marker heatmap: per-state mean, z per gene across states
mean = E.groupby("LS")[present].mean()
z = (mean - mean.mean(0)) / (mean.std(0) + 1e-9)
z.to_csv(os.path.join(D2, "state_annotation_heatmap_z.csv"))

# per-state cross-cytogenetic breadth (effective # subtypes) + patient dominance
ann = pd.read_csv(C.ANNOT).set_index("LS")
rows = []
for s in z.index:
    cells = lscol[lscol == s].index
    cv = cyto.loc[cells].value_counts(normalize=True); cv = cv[cv > 0]
    eff = float(np.exp(-(cv * np.log(cv)).sum()))          # Hill q=1 (effective # subtypes)
    npat = int(samp.loc[cells].nunique())
    toppt = float(samp.loc[cells].value_counts(normalize=True).iloc[0])
    cls = "single-pt" if toppt >= 0.77 else ("few-pt" if npat <= 3 else "multi-pt")
    rows.append([s, ann.lineage_program.get(s, "unresolved"), int((lscol == s).sum()), round(eff, 2), npat, round(toppt, 2), cls])
meta = pd.DataFrame(rows, columns=["LS","lineage_program","n_cells","eff_n_cyto","n_patients","top_pt_frac","patient_class"])
meta.to_csv(os.path.join(D2, "state_annotation_meta.csv"), index=False)
print(f"wrote state_annotation_heatmap_z.csv ({z.shape[0]} states x {z.shape[1]} genes) and meta")
print(meta.patient_class.value_counts().to_string())
print("proliferative states:", meta.LS[meta.lineage_program == "proliferative"].tolist())
