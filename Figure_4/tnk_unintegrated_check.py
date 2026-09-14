#!/usr/bin/env python
"""Do AML T/NK cells co-cluster with healthy-BM T/NK cells, or separate? Re-embed the SAME T/NK cells
WITHOUT batch correction (un-integrated PCA) and quantify overlap. If AML and HBM still mix -> AML T/NK
are transcriptionally normal-like; if they separate, also report per-cluster patient count to tell
AML-reprogramming from mere per-patient batch. No fabricated values."""
import os, numpy as np, pandas as pd, scanpy as sc
H5="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
D ="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
TNK=["AML-CD4T","AML-CTL","AML-NK","AML-Naïve_CD8T","Activated_CD4T","CTL","NK","Naïve_CD4T","Naïve_CD8T"]
a=sc.read_h5ad(H5); a=a[a.obs["Cell_Type"].astype(str).isin(TNK)].copy()
a.X=a.layers["Raw_Counts"].copy(); sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a); a.raw=a
sc.pp.highly_variable_genes(a,n_top_genes=2000); a=a[:,a.var.highly_variable].copy()
sc.pp.scale(a,max_value=10); sc.tl.pca(a,n_comps=30)
sc.pp.neighbors(a,n_neighbors=15,use_rep="X_pca")               # NO harmony
sc.tl.leiden(a,resolution=0.8,key_added="cl",flavor="igraph",n_iterations=2,directed=False)
sc.tl.umap(a)
grp=np.where(a.obs["SampleID"].astype(str).str.startswith("0_HealthyBM"),"HBM","AML")
out=pd.DataFrame({"barcode":a.obs_names,"SampleID":a.obs["SampleID"].astype(str).values,
                  "group":grp,"cl":a.obs["cl"].values,
                  "UMAP1":a.obsm["X_umap"][:,0],"UMAP2":a.obsm["X_umap"][:,1]})
out.to_csv(os.path.join(D,"tnk_percell_unintegrated.csv"),index=False)
# per-cluster composition + patient diversity
t=out.groupby("cl")["group"].value_counts().unstack(fill_value=0)
t["HBM_frac"]=t.get("HBM",0)/t.sum(axis=1)
t["n_AML_pts"]=out[out.group=="AML"].groupby("cl")["SampleID"].nunique().reindex(t.index).fillna(0).astype(int)
t["top_AML_pt_share"]=out[out.group=="AML"].groupby("cl")["SampleID"].apply(
    lambda s: s.value_counts(normalize=True).iloc[0] if len(s) else np.nan).reindex(t.index).round(2)
print("un-integrated T/NK: %d cells, %d clusters | HBM overall frac %.3f"%(len(out),out.cl.nunique(),(grp=="HBM").mean()))
print("\nper-cluster (un-integrated): AML/HBM counts, HBM_frac, #AML patients, largest single-patient share of the AML cells")
print(t.sort_index(key=lambda x:x.astype(int)).to_string())
mixed=((t["HBM_frac"]>0.1)&(t["HBM_frac"]<0.9)).mean()
print("\nfraction of clusters that are MIXED (10-90%% HBM): %.2f"%mixed)
print("clusters ~pure AML (HBM<5%%):",[i for i in t.index if t.loc[i,"HBM_frac"]<0.05])
print("clusters ~pure HBM (HBM>95%%):",[i for i in t.index if t.loc[i,"HBM_frac"]>0.95])
