#!/usr/bin/env python
"""Re-embed and sub-cluster the LS10 leukemic cells (AML samples) to probe internal heterogeneity.
Fresh normalize -> HVG -> PCA -> neighbors -> UMAP -> Leiden. Diagnose whether substructure is
cell-cycle phase, patient/batch, or a real program (crosstabs + marker/covariate UMAPs). Saves
subcluster labels + UMAP for downstream. No fabricated values."""
import os, sys, numpy as np, pandas as pd, scanpy as sc, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"); import config as C
OUT = os.path.join(C.BASE, "04_Main_Figures/Figure_3_panels"); D3 = C.DATA
a = sc.read_h5ad(C.H5, backed="r")
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
lsser = pd.Series([bc2ls.get(b) for b in a.obs_names], index=a.obs_names)
mask = (lsser.values == "LS_10") & (a.obs["SampleType"].astype(str) == "AML").values
s = a[mask].to_memory()
print("LS10 AML cells:", s.n_obs)
s.X = s.layers["Raw_Counts"].copy()
sc.pp.normalize_total(s, target_sum=1e4); sc.pp.log1p(s)
s.raw = s
sc.pp.highly_variable_genes(s, n_top_genes=2000); s = s[:, s.var.highly_variable]
sc.pp.scale(s, max_value=10)
sc.tl.pca(s, n_comps=30)
sc.pp.neighbors(s, n_neighbors=15, n_pcs=30)
sc.tl.umap(s)
sc.tl.leiden(s, resolution=0.5, key_added="ls10_sub")
print("subclusters:", s.obs["ls10_sub"].value_counts().sort_index().to_dict())
# composition diagnostics
print("\nsubcluster x phase (row %):")
print(pd.crosstab(s.obs["ls10_sub"], s.obs["phase"], normalize="index").round(2).to_string())
print("\nsubcluster x patient (counts):")
print(pd.crosstab(s.obs["ls10_sub"], s.obs["SampleID"]).to_string())
# save labels + umap
out = pd.DataFrame({"barcode": s.obs_names, "ls10_sub": s.obs["ls10_sub"].values,
                    "phase": s.obs["phase"].values, "SampleID": s.obs["SampleID"].astype(str).values,
                    "UMAP1": s.obsm["X_umap"][:,0], "UMAP2": s.obsm["X_umap"][:,1]})
out.to_csv(os.path.join(D3, "LS10_subcluster.csv"), index=False)
# figure: UMAP colored by subcluster / phase / patient / key markers
for g in ["MKI67","MPO","CD34","S_score","G2M_score"]:
    if g in ("S_score","G2M_score"): continue
sc.settings.figdir = OUT
fig = sc.pl.umap(s, color=["ls10_sub","phase","SampleID","MKI67","MPO","CD34"], ncols=3,
                 wspace=0.35, show=False, return_fig=True)
fig.savefig(os.path.join(OUT, "Figure_3_LS10_subcluster_umap.png"), dpi=160, bbox_inches="tight")
print("wrote LS10_subcluster.csv and Figure_3_LS10_subcluster_umap.png")
