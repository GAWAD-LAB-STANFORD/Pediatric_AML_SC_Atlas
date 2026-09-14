#!/usr/bin/env python
"""Figure 3C data: per-cell Wilcoxon DE of LS10 (poor cell-cycle) vs the non-prognostic
proliferation comparators (LS39, LS3). Writes the full DE plus the non-histone up-genes
(input to the GO panel) and the down-genes (up in comparators). Histones are split out so
the GO panel reflects program, not S-phase replication-histone load. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py, scanpy as sc, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
D3 = C.DATA
poor = ["LS_10"]; comp = ["LS_39", "LS_3"]
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
ra["LS"] = ra["res4.50"].map(raw2ls)
sub = ra[ra.LS.isin(poor + comp)].copy(); sub["grp"] = np.where(sub.LS.isin(poor), "LS10_poor", "comp_ns")
bc2grp = dict(zip(sub.barcode, sub.grp))
f = h5py.File(C.H5, "r"); names = [x.decode() if isinstance(x, bytes) else x for x in f["obs"]["_index"][:]]; f.close()
idx = np.where(np.isin(names, sub.barcode.values))[0]
ad = sc.read_h5ad(C.H5, backed="r")[idx].to_memory()
ad.X = ad.layers["Raw_Counts"].copy(); sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
ad.obs["grp"] = [bc2grp.get(b, "NA") for b in ad.obs_names]; ad = ad[ad.obs.grp != "NA"].copy()
print("cells:", ad.obs.grp.value_counts().to_dict())
sc.tl.rank_genes_groups(ad, "grp", groups=["LS10_poor"], reference="comp_ns", method="wilcoxon")
r = sc.get.rank_genes_groups_df(ad, group="LS10_poor")
r.to_csv(os.path.join(D3, "LS10_vs_prolifNS_DE.csv"), index=False)
sig = r[r.pvals_adj < 0.05]
up = sig[sig.logfoldchanges > 0.5].copy(); up["histone"] = up.names.str.match(r"^HIST\d")
dn = sig[sig.logfoldchanges < -0.5].sort_values("logfoldchanges")
up[~up.histone][["names", "logfoldchanges", "pvals_adj"]].to_csv(os.path.join(D3, "LS10_up_nonhistone.csv"), index=False)
dn[["names", "logfoldchanges", "pvals_adj"]].to_csv(os.path.join(D3, "LS10_dn.csv"), index=False)
print(f"UP {len(up)} ({int(up.histone.sum())} histone) -> {int((~up.histone).sum())} non-histone; DOWN {len(dn)}")
print("wrote LS10_vs_prolifNS_DE.csv, LS10_up_nonhistone.csv, LS10_dn.csv")
