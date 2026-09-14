#!/usr/bin/env python
"""Per-column (49 leukemic states + 8 normal compartments) module scores from log-normalized
expression (h5ad 'counts' layer = log1p CP10k). Two products:
  (1) proliferation module score  -> top annotation strip for the Fig 2 regulon panel
  (2) canonical lineage-marker module scores -> evidence to label the 'unresolved' states
Columns are defined IDENTICALLY to compute_regulon_landscape.py (leukemic via REFA res4.50 ->
stable LS; normal via Cell_Type). Mean per column, then z across columns. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
D3 = C.DATA

PROLIF = ["MKI67","TOP2A","PCNA","CCNB1","CCNB2","CDK1","MCM2","MCM6","TYMS",
          "BIRC5","CENPF","UBE2C","AURKB","CCNA2","RRM2","CDC20","CENPA","CKS1B"]
LIN = {
 "HSC/MPP":     ["CRHBP","HLF","AVP","MEIS1","HOXA9","MLLT3","GATA2","ERG","MECOM","PRSS57"],
 "GMP":         ["MPO","ELANE","PRTN3","AZU1","CTSG","RNASE2","CEBPE","CEBPA"],
 "Monocyte":    ["CD14","FCN1","VCAN","S100A8","S100A9","LYZ","CSF1R","CD68","MNDA"],
 "cDC":         ["FCER1A","CD1C","CLEC10A","CLEC9A","XCR1","ITGAX"],
 "pDC":         ["IL3RA","LILRA4","IRF7","GZMB","TCF4","CLEC4C"],
 "Erythroid":   ["GATA1","KLF1","HBB","HBA1","ALAS2","GYPA","TFRC","CA1"],
 "Megakaryocyte":["PF4","PPBP","ITGA2B","GP9","VWF","PLEK","GP1BA"],
 "B/plasma":    ["CD19","MS4A1","CD79A","CD79B","EBF1","VPREB1","MZB1","IGHM","IGLL1"],
 "T/NK":        ["CD3D","CD3E","CD2","IL7R","NKG7","GNLY","KLRD1","CCL5","TCF7"],
 "Neutrophil":  ["FUT4","CSF3R","FCGR3B","CEACAM8","MMP9","S100A12","CXCR2"],
}
NORM = {"0_HSPC":"NORM_HSPC","Myeloid_Pro":"NORM_MyeloidPro",
        "CD14_Monocyte":"NORM_Monocyte","CD16_Monocyte":"NORM_Monocyte","Macrophage":"NORM_Macrophage",
        "mDC":"NORM_DC","pDC":"NORM_DC","Erythrocytes":"NORM_Erythroid",
        "Naïve_CD4T":"NORM_T_NK","Naïve_CD8T":"NORM_T_NK","CTL":"NORM_T_NK","NK":"NORM_T_NK","Activated_CD4T":"NORM_T_NK",
        "CD20+B":"NORM_B","ProB":"NORM_B","PreB":"NORM_B","PlasmaB":"NORM_B","CD34+ProB":"NORM_B"}

f = h5py.File(C.H5, "r")
var_names = np.array([v.decode() if isinstance(v, bytes) else v
                      for v in (f["var"]["_index"][:] if "_index" in f["var"] else f["var"][list(f["var"].keys())[0]][:])])
obs_names = np.array([o.decode() if isinstance(o, bytes) else o
                      for o in (f["obs"]["_index"][:] if "_index" in f["obs"] else f["obs"][list(f["obs"].keys())[0]][:])])
ct_codes = f["obs"]["Cell_Type"]["codes"][:]; ct_cats = np.array([c.decode() if isinstance(c, bytes) else c for c in f["obs"]["Cell_Type"]["categories"][:]])
cell_type = pd.Series(ct_cats[ct_codes], index=obs_names)
allg = [g for g in (PROLIF + [x for v in LIN.values() for x in v]) if g in set(var_names)]
allg = list(dict.fromkeys(allg)); gpos = {g: i for i, g in enumerate(var_names)}
gi = [gpos[g] for g in allg]
grp = f["layers"]["counts"]; shape = tuple(grp.attrs["shape"])
Xfull = csr_matrix((grp["data"][:], grp["indices"][:], grp["indptr"][:]), shape=shape)  # ~1.6 GB
f.close()
Xsub = np.asarray(Xfull[:, gi].todense())   # cells x marker genes (log-normalized)
E = pd.DataFrame(Xsub, index=obs_names, columns=allg)

# ---- column assignment: leukemic (REFA -> LS) then normal (Cell_Type -> NORM) ----
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
ctd = cell_type.to_dict()
col = []
for b in E.index:
    v = bc2ls.get(b)
    if isinstance(v, str) and v.startswith("LS_"): col.append(v)
    else: col.append(NORM.get(ctd.get(b, ""), None))
E["col"] = col; E = E[E["col"].notna()]

colmean = E.groupby("col")[allg].mean()   # per-column per-gene mean log-norm expression
colmean.to_csv(os.path.join(D3, "marker_gene_colmean.csv"))
def modscore(genes):
    g = [x for x in genes if x in E.columns]
    return E.groupby("col")[g].mean().mean(1)   # mean log-norm expr across module genes, per column

# proliferation strip (mean, and z across columns)
prolif_mean = modscore(PROLIF)
prolif_z = (prolif_mean - prolif_mean.mean()) / (prolif_mean.std() + 1e-9)
pd.DataFrame({"proliferation_mean": prolif_mean, "proliferation_z": prolif_z}).to_csv(
    os.path.join(D3, "proliferation_strip.csv"))

# lineage module matrix (z across columns) -> assignment evidence
L = pd.DataFrame({k: modscore(v) for k, v in LIN.items()})
Lz = (L - L.mean(0)) / (L.std(0) + 1e-9)
Lz.to_csv(os.path.join(D3, "lineage_module_z.csv"))

ann = pd.read_csv(C.ANNOT).set_index("LS")
rows = []
for s in [c for c in Lz.index if c.startswith("LS_")]:
    srt = Lz.loc[s].sort_values(ascending=False)
    rows.append([s, ann.lineage_program.get(s, ""), round(prolif_z[s], 2),
                 srt.index[0], round(srt.iloc[0], 2), srt.index[1], round(srt.iloc[1], 2),
                 round(srt.iloc[0] - srt.iloc[1], 2)])
d = pd.DataFrame(rows, columns=["LS","current","prolif_z","top_lineage","z1","2nd","z2","margin"])
d.to_csv(os.path.join(D3, "lineage_marker_assignment.csv"), index=False)
pd.set_option("display.width", 170); pd.set_option("display.max_rows", 60)
print("== UNRESOLVED states: marker-module evidence (z across columns) ==")
print(d[d.current == "unresolved"].to_string(index=False))
print("\n== proliferative states ==")
print(d[d.current == "proliferative"].to_string(index=False))
print("\nproliferation strip (top 8 by z):")
print(prolif_z.sort_values(ascending=False).head(8).round(2).to_string())
print("wrote proliferation_strip.csv, lineage_module_z.csv, lineage_marker_assignment.csv")
