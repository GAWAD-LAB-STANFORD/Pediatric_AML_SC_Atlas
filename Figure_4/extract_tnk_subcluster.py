#!/usr/bin/env python
"""Figure 4 rebuild, step 1: re-derive finer T/NK-cell states from the ORIGINAL count matrix.
The current h5ad only stores coarse T/NK labels (CTL, NK, CD4T, Naive_CD8T); the finer subsets the
figure needs (GZMB/GZMK-CD8T, MAIT, effector-memory, Treg, etc.) are not stored, so we pull the T/NK
cells out of the raw counts and re-subcluster them. Harmony-integrated by SampleID so shared T/NK states
merge across patients. Writes per-cell subcluster labels + UMAP + per-cluster marker table for annotation.
No fabricated values."""
import os, numpy as np, pandas as pd, scanpy as sc, anndata as ad
H5  = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
OUT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
os.makedirs(OUT, exist_ok=True)
sc.settings.verbosity = 1
TNK = ["AML-CD4T","AML-CTL","AML-NK","AML-Naïve_CD8T","Activated_CD4T","CTL","NK","Naïve_CD4T","Naïve_CD8T"]
MARKERS = ["CD3D","CD3E","CD8A","CD8B","CD4","IL7R","CCR7","SELL","TCF7",   # T backbone / naive
           "GZMK","GZMB","GZMA","GNLY","PRF1","NKG7","KLRD1","NCAM1","FCGR3A",  # cytotoxic / NK
           "SLC4A10","KLRB1","TRAV1-2","FOXP3","IL2RA","CTLA4","MKI67","CD14","LYZ"]  # MAIT / Treg / prolif / myeloid-contam

adata = sc.read_h5ad(H5)
print("full:", adata.shape)
ct = adata.obs["Cell_Type"].astype(str)
adata = adata[ct.isin(TNK)].copy()
print("T/NK subset:", adata.shape)
# rebuild from raw counts
adata.X = adata.layers["Raw_Counts"].copy()
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
adata.raw = adata
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=30)
# batch-integrate by patient so shared T/NK states merge
try:
    sc.external.pp.harmony_integrate(adata_hvg, "SampleID")
    rep = "X_pca_harmony"
except Exception as e:
    print("harmony failed, using PCA:", e); rep = "X_pca"
sc.pp.neighbors(adata_hvg, n_neighbors=15, use_rep=rep)
sc.tl.leiden(adata_hvg, resolution=0.8, key_added="tnk_leiden", flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(adata_hvg)
adata.obs["tnk_leiden"] = adata_hvg.obs["tnk_leiden"].values
adata.obsm["X_umap_tnk"] = adata_hvg.obsm["X_umap"]
print("clusters:", adata.obs["tnk_leiden"].value_counts().sort_index().to_dict())
# per-cluster marker means (log-norm) for annotation
present = [g for g in MARKERS if g in adata.raw.var_names]
E = pd.DataFrame(adata.raw[:, present].X.toarray(), columns=present, index=adata.obs_names)
E["cl"] = adata.obs["tnk_leiden"].values
mean_by_cl = E.groupby("cl")[present].mean()
mean_by_cl.to_csv(os.path.join(OUT, "tnk_cluster_marker_means.csv"))
# save per-cell table (barcode, sample, coarse type, subcluster, UMAP)
out = pd.DataFrame({"barcode": adata.obs_names,
                    "SampleID": adata.obs["SampleID"].astype(str).values,
                    "coarse_type": ct.reindex(adata.obs_names).values if False else adata.obs["Cell_Type"].astype(str).values,
                    "tnk_leiden": adata.obs["tnk_leiden"].values,
                    "UMAP1": adata.obsm["X_umap_tnk"][:,0], "UMAP2": adata.obsm["X_umap_tnk"][:,1]})
out.to_csv(os.path.join(OUT, "tnk_percell.csv"), index=False)
print("\nper-cluster marker means (for annotation):")
print(mean_by_cl.round(2).to_string())
print("\nwrote tnk_percell.csv + tnk_cluster_marker_means.csv to", OUT)
