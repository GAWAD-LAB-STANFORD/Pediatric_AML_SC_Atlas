#!/usr/bin/env python
"""LS19 characterization data for the pivoted Figure 3: per-group mean log-norm expression + %
expressing of pDC, progenitor/stem, interferon, differentiation, and surface-target markers, across
LS19 vs normal pDC vs normal HSPC vs other leukemic. Tests the 'least-differentiated pDC-committed'
story. CSR-direct (fast). Also LS19 TARGET fraction by cytogenetic subtype. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"); import config as C
D3 = C.DATA
GENES = {
 "pDC":         ["IRF7","LILRA4","GZMB","TCF4","CLEC4C","IRF8","SPIB","BST2"],
 "progenitor":  ["SPINK2","PRSS57","FAM30A","SPNS3","CD34","HOXA9","MEIS1","CD38"],
 "interferon":  ["IFI6","IFITM1","XAF1","ISG15","MX1"],
 "differentiation": ["MPO","ELANE","LYZ","CD14","CD1C","ITGAX"],
 "surface":     ["IL3RA","CLEC4C","BST2","CD96","CD33","FLT3"],
}
mk = list(dict.fromkeys([g for v in GENES.values() for g in v]))
f = h5py.File(C.H5, "r"); dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
obs = dec(f["obs"]["_index"][:]); var = dec(f["var"]["_index"][:])
ct = pd.Series(dec(f["obs"]["Cell_Type"]["categories"][:])[f["obs"]["Cell_Type"]["codes"][:]], index=obs)
st = pd.Series(dec(f["obs"]["SampleType"]["categories"][:])[f["obs"]["SampleType"]["codes"][:]], index=obs)
gpos = {g: i for i, g in enumerate(var)}; present = [g for g in mk if g in gpos]; gi = [gpos[g] for g in present]
grp = f["layers"]["counts"]; X = csr_matrix((grp["data"][:], grp["indices"][:], grp["indptr"][:]), shape=tuple(grp.attrs["shape"]))
f.close()
E = pd.DataFrame(np.asarray(X[:, gi].todense()), columns=present, index=obs)
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
lsser = pd.Series([bc2ls.get(b) for b in obs], index=obs)
def group(b):
    v = bc2ls.get(b)
    if v == "LS_19": return "LS19 (poor)"
    if isinstance(v, str) and v.startswith("LS_"): return "other leukemic"
    if ct.get(b) == "pDC": return "normal pDC"
    if ct.get(b) == "0_HSPC": return "normal HSPC"
    return None
E["grp"] = [group(b) for b in obs]; E = E[E["grp"].notna()]
order = ["normal HSPC", "LS19 (poor)", "normal pDC", "other leukemic"]
mean = E.groupby("grp")[present].mean().reindex(order)
pct = E.groupby("grp")[present].apply(lambda d: (d > 0).mean()).reindex(order)
# long format for dotplot
rows = []
for prog, gg in GENES.items():
    for g in gg:
        if g not in present: continue
        for grp_ in order:
            rows.append([prog, g, grp_, mean.loc[grp_, g], pct.loc[grp_, g]])
out = pd.DataFrame(rows, columns=["program","gene","group","mean","pct"])
out.to_csv(os.path.join(D3, "LS19_marker_dotplot.csv"), index=False)
print("group sizes:", E["grp"].value_counts().to_dict())
print("\nkey markers (mean log-norm) — LS19 vs refs:")
for g in ["SPINK2","PRSS57","CD34","IL3RA","CLEC4C","BST2","IRF7","LILRA4","IFI6","MPO","LYZ"]:
    if g in present: print(f"  {g:8} HSPC={mean.loc['normal HSPC',g]:.2f} LS19={mean.loc['LS19 (poor)',g]:.2f} pDC={mean.loc['normal pDC',g]:.2f} otherLeuk={mean.loc['other leukemic',g]:.2f}")
