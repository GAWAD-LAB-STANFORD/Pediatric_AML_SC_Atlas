#!/usr/bin/env python
"""Figure 3F data: LS10 surface-marker screen. (1) DE of LS10 vs all other leukemic cells
restricted to the Bausch-Fluck in-silico surfaceome. (2) HSC-sparing candidates = surface
genes present on LS10 (>=30%) but absent on gated normal HSC (<=5%) and full HSPC (<=10%).
(3) Dot-plot data (% + mean) for the top HSC-sparing candidates across cell groups.
Reference-subsampled for the DE. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py, scanpy as sc, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
D3 = C.DATA; D2 = D3.replace("Figure_3", "Figure_2")
SURF = set(l.strip() for l in open(os.path.join(C.BASE, "ref/surfaceome_genes.txt")) if l.strip())
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
ra["LS"] = ra["res4.50"].map(raw2ls)
gated = set(pd.read_csv(os.path.join(D2, "normal_HSC_gated_barcodes.csv"))["barcode"].astype(str))
rng = np.random.default_rng(0)
f = h5py.File(C.H5, "r")
names = np.array([x.decode() if isinstance(x, bytes) else x for x in f["obs"]["_index"][:]])
nd = f["obs"]["Cell Type"]; cc = [x.decode() if isinstance(x, bytes) else x for x in nd["categories"][:]]
ctype = np.array([cc[i] if i >= 0 else "NA" for i in nd["codes"][:]], dtype=object); f.close()

# --- (1)+(2) LS10 vs other leukemic (subsample reference), surfaceome ---
oth = ra.barcode[ra.LS != "LS_10"].values; oth_s = set(rng.choice(oth, 8000, replace=False))
keep = set(ra.barcode[ra.LS == "LS_10"]) | oth_s | gated | set(names[ctype == "0_HSPC"])
bc2 = {}
for b, l in zip(ra.barcode, ra.LS):
    if l == "LS_10": bc2[b] = "LS10"
    elif b in oth_s: bc2[b] = "other_leuk"
for b in gated: bc2[b] = "HSC_gated"
idx = np.where(np.isin(names, list(keep)))[0]
name2ct = dict(zip(names, ctype))
ad = sc.read_h5ad(C.H5, backed="r")[idx].to_memory()
ad.X = ad.layers["Raw_Counts"].copy(); sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
ad.obs["grp"] = [bc2.get(b, "HSPC_full" if name2ct.get(b) == "0_HSPC" else "NA") for b in ad.obs_names]
de = ad[ad.obs.grp.isin(["LS10", "other_leuk"])].copy()
sc.tl.rank_genes_groups(de, "grp", groups=["LS10"], reference="other_leuk", method="wilcoxon")
r = sc.get.rank_genes_groups_df(de, group="LS10")
r = r[(r.names.isin(SURF)) & (r.pvals_adj < 0.05) & (r.logfoldchanges > 0.5)].copy()
r.sort_values("logfoldchanges", ascending=False).to_csv(os.path.join(D3, "LS10_surface_markers.csv"), index=False)

def pct(mask, gene):
    x = ad[:, gene].X; x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    return (x[mask] > 0).mean() * 100
sg = [g for g in r.names if g in ad.var_names]
p = pd.DataFrame({"gene": sg,
    "pct_LS10": [pct((ad.obs.grp == "LS10").values, g) for g in sg],
    "pct_HSC_gated": [pct((ad.obs.grp == "HSC_gated").values, g) for g in sg],
    "pct_HSPC_full": [pct(np.isin(ad.obs.grp, ["HSC_gated", "HSPC_full"]), g) for g in sg]})
cand = p[(p.pct_LS10 >= 30) & (p.pct_HSC_gated <= 5) & (p.pct_HSPC_full <= 10)].sort_values("pct_LS10", ascending=False)
cand.to_csv(os.path.join(D3, "LS10_surface_vs_HSPC_candidates.csv"), index=False)
print("HSC-sparing surface candidates:\n", cand.round(1).to_string(index=False))

# --- (3) dot-plot data for the candidate markers across groups (incl. LS39/LS3) ---
markers = ["CD96", "CLEC12A", "CD68", "SUCNR1"]
bc3 = {}
for b, l in zip(ra.barcode, ra.LS):
    if l == "LS_10": bc3[b] = "LS10 (poor)"
    elif l == "LS_39": bc3[b] = "LS39 (ns)"
    elif l == "LS_3": bc3[b] = "LS3 (ns)"
    elif b in oth_s: bc3[b] = "other leukemic"
for b in gated: bc3[b] = "normal HSC"
ad.obs["g3"] = [bc3.get(b, "NA") for b in ad.obs_names]
rows = []
for m in markers:
    if m not in ad.var_names: continue
    x = ad[:, m].X; x = np.asarray(x.todense()).ravel() if hasattr(x, "todense") else np.asarray(x).ravel()
    for g in ["normal HSC", "other leukemic", "LS3 (ns)", "LS39 (ns)", "LS10 (poor)"]:
        mm = (ad.obs.g3 == g).values
        if mm.sum(): rows.append([m, g, round((x[mm] > 0).mean() * 100, 1), round(x[mm].mean(), 3), int(mm.sum())])
pd.DataFrame(rows, columns=["gene", "group", "pct", "mean_expr", "n"]).to_csv(os.path.join(D3, "LS10_surface_dotplot.csv"), index=False)
print("wrote LS10_surface_markers.csv, LS10_surface_vs_HSPC_candidates.csv, LS10_surface_dotplot.csv")
