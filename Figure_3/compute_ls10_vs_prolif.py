#!/usr/bin/env python
"""Figure 3C data: what distinguishes LS10 among proliferating blasts. Differential expression of
LS10 vs the other three proliferation-high states (LS3, LS20, LS39), restricted to AML-sample
LEUKEMIC cells only (healthy-BM cells that carry a proliferation Cell_Type label and cluster into
these states are excluded, so the contrast is purely leukemic). Wilcoxon, BH-adjusted; per-group
percent expressing. Writes LS10_vs_prolif_DE.csv. No fabricated values."""
import sys, numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"); import config as C
a = sc.read_h5ad(C.H5, backed="r")
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
lsser = pd.Series([bc2ls.get(b) for b in a.obs_names], index=a.obs_names)
prolif = ["LS_3", "LS_10", "LS_20", "LS_39"]
is_aml = (a.obs["SampleType"].astype(str) == "AML").values          # exclude healthy-BM cells
mask = lsser.isin(prolif).values & is_aml
n_excl = int((lsser.isin(prolif).values & ~is_aml).sum())
sub = a[mask].to_memory(); sub.obs["LS"] = lsser[mask].values; sub.X = sub.layers["counts"]
sub.obs["grp"] = np.where(sub.obs["LS"] == "LS_10", "LS10", "other_prolif")
print("cells (AML-only):", sub.obs["LS"].value_counts().to_dict(), "| healthy-BM cells excluded:", n_excl)
sc.tl.rank_genes_groups(sub, "grp", groups=["LS10"], reference="other_prolif", method="wilcoxon")
r = sc.get.rank_genes_groups_df(sub, group="LS10")
X = sub.X; X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
is10 = (sub.obs["grp"] == "LS10").values; gi = {g: i for i, g in enumerate(sub.var_names)}
r["pct_LS10"]  = [float((X[is10, gi[g]] > 0).mean())  for g in r.names]
r["pct_other"] = [float((X[~is10, gi[g]] > 0).mean()) for g in r.names]
r.to_csv(C.DATA + "/LS10_vs_prolif_DE.csv", index=False)
up = r[(r.pvals_adj<0.05)&(r.logfoldchanges>0.5)&(r.pct_LS10>0.25)].sort_values("logfoldchanges",ascending=False).head(10)
dn = r[(r.pvals_adj<0.05)&(r.logfoldchanges<-0.5)&(r.pct_other>0.25)].sort_values("logfoldchanges").head(10)
print("top10 UP:", list(up.names)); print("top10 DOWN:", list(dn.names))
print("wrote LS10_vs_prolif_DE.csv (AML-sample leukemic cells only)")
