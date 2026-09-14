#!/usr/bin/env python
"""Supplementary-figure data: strict true-LSC gating and characterization.
LSC definition (literature: Ding 2017; Kikushige PNAS 2010 for TIM3; CLEC12A/CD96) =
CD34+ CD38- AND quiescent (HLF/AVP+) AND non-cycling (S<0 & G2M<0) AND carrying an
LSC-aberrant marker (TIM3/CD96/CLEC12A) that normal HSC lack. Outputs: the gating funnel
(AML vs healthy BM), the strict-LSC barcodes (for the UMAP), and the marker profile of the
candidates vs normal HSC (primitive) and normal GMP (Myeloid_Pro, committed). No fabrication."""
import os, sys, numpy as np, pandas as pd, h5py, scanpy as sc, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
D3 = C.DATA; D2 = D3.replace("Figure_3", "Figure_2")
f = h5py.File(C.H5, "r")
def cat(n):
    nd = f["obs"][n]; c = [x.decode() if isinstance(x, bytes) else x for x in nd["categories"][:]]
    return np.array([c[i] if i >= 0 else "NA" for i in nd["codes"][:]], dtype=object)
names = np.array([x.decode() if isinstance(x, bytes) else x for x in f["obs"]["_index"][:]])
stype = cat("SampleType"); ctype = cat("Cell Type"); sid = cat("SampleID")
S = f["obs"]["S_score"][:]; G2M = f["obs"]["G2M_score"][:]
genes = np.array([g.decode() if isinstance(g, bytes) else g for g in f["var"]["_index"][:]])
gi = {g: int(np.where(genes == g)[0][0]) for g in ["HLF", "AVP", "CD34", "CD38", "HAVCR2", "CD96", "CLEC12A"]}
RC = f["layers"]["Raw_Counts"]; ind = RC["indices"]; ip = RC["indptr"][:]; N = len(names)
def present(g):
    tgt = gi[g]; hits = []
    for a in range(0, ind.shape[0], 30_000_000):
        ch = ind[a:a+30_000_000]; loc = np.where(ch == tgt)[0]
        if len(loc): hits.append(loc + a)
    p = np.concatenate(hits) if hits else np.array([], int)
    m = np.zeros(N, bool); m[np.searchsorted(ip, p, side="right") - 1] = True; return m
P = {g: present(g) for g in gi}; f.close()
cd34p, cd38n = P["CD34"], ~P["CD38"]; quies = P["HLF"] | P["AVP"]; noncyc = (S < 0) & (G2M < 0)
lsc_ab = P["HAVCR2"] | P["CD96"] | P["CLEC12A"]
# funnel
rows = []
for lab, grp in [("AML sample", stype == "AML"), ("healthy BM", stype == "HealthyBM")]:
    a = grp & cd34p; b = a & cd38n; c = b & quies; d = c & noncyc; e = d & lsc_ab
    for step, mask in [("CD34+", a), ("+ CD38-", b), ("+ quiescent (HLF/AVP+)", c), ("+ non-cycling", d), ("+ LSC marker (TIM3/CD96/CLEC12A)", e)]:
        rows.append([lab, step, int(mask.sum())])
pd.DataFrame(rows, columns=["group", "step", "count"]).to_csv(os.path.join(D3, "strict_LSC_funnel.csv"), index=False)
strict = (stype == "AML") & cd34p & cd38n & quies & noncyc & lsc_ab
pd.DataFrame({"barcode": names[strict]}).to_csv(os.path.join(D3, "strict_LSC_barcodes.csv"), index=False)
print(f"strict true-LSC: AML {int(strict.sum())} vs healthyBM {int(((stype=='HealthyBM')&cd34p&cd38n&quies&noncyc&lsc_ab).sum())}")
print("  patients:", dict(pd.Series(sid[strict]).value_counts()))
# marker profile: candidates vs gated HSC vs normal GMP (Myeloid_Pro)
gated = set(pd.read_csv(os.path.join(D2, "normal_HSC_gated_barcodes.csv"))["barcode"].astype(str))
grp = np.array([None] * N, dtype=object)
grp[np.array([b in gated for b in names])] = "normal HSC (primitive)"
grp[ctype == "Myeloid_Pro"] = "normal GMP (committed)"
grp[strict] = "strict-LSC candidates (n=%d)" % int(strict.sum())
sel = np.where(grp != None)[0]
ad = sc.read_h5ad(C.H5, backed="r")[sel].to_memory()
ad.X = ad.layers["Raw_Counts"].copy(); sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
ad.obs["grp"] = grp[sel]
MOD = {"stem/quiescence": ["CD34", "HLF", "AVP", "HOXA9"], "LSC-aberrant": ["CD96", "CLEC12A", "HAVCR2"],
       "GMP/granule": ["MPO", "ELANE", "AZU1", "CTSG"], "monocyte": ["LYZ", "CD14"]}
out = []
for mod, gs in MOD.items():
    for g in gs:
        if g in ad.var_names:
            x = ad[:, g].X; x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
            for gr in ad.obs.grp.unique():
                m = (ad.obs.grp == gr).values
                out.append([g, mod, gr, round((x[m] > 0).mean() * 100, 1), round(x[m].mean(), 3)])
pd.DataFrame(out, columns=["gene", "module", "group", "pct", "mean_expr"]).to_csv(os.path.join(D3, "strictLSC_marker_profile.csv"), index=False)
print("wrote strict_LSC_funnel.csv, strict_LSC_barcodes.csv, strictLSC_marker_profile.csv")
