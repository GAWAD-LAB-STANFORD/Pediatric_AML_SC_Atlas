#!/usr/bin/env python
"""Figure 3 step 4: per-cell Wilcoxon DE for the leukemic states (core cells), each
state vs the rest -> up-gene lists (padj<0.001 & log2FC>1) for GO:BP (panel E)."""
import os, sys
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"): os.environ[v]="4"
import scanpy as sc, numpy as np, pandas as pd, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
f=h5py.File(C.H5,"r"); nd=f["obs"]["Cell Type"]
cats=[c.decode() if isinstance(c,bytes) else c for c in nd["categories"][:]]
ct=np.array([cats[i] if i>=0 else "NA" for i in nd["codes"][:]],dtype=object); f.close()
mask=np.isin(ct,C.LEUK_CELLTYPES); pos=np.where(mask)[0]
ref=pd.read_csv(C.REFA)["res4.50"].astype(int).values
raws,ls=C.stable_raws(); raw2ls={rw:i for i,rw in enumerate(raws)}
is_core=np.isin(ref,raws); keep=pos[is_core]; lslab=[ls[raw2ls[r]] for r in ref[is_core]]
print(f"core leukemic cells={len(keep)}, states={len(ls)}",flush=True)
ad=sc.read_h5ad(C.H5,backed="r")[keep].to_memory(); ad.obs["LS"]=pd.Categorical(lslab,categories=ls)
mn=float(ad.X.min()); mx=float(ad.X.max())
if mn<-0.01:
    lyr="counts" if "counts" in ad.layers else "Raw_Counts"; ad.X=ad.layers[lyr].copy()
    sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
elif mx>50:
    sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad)
sc.tl.rank_genes_groups(ad,"LS",method="wilcoxon")
r=ad.uns["rank_genes_groups"]; out=[]
for g in r["names"].dtype.names:
    out.append(pd.DataFrame({"LS":g,"gene":r["names"][g],"log2FC":r["logfoldchanges"][g],"padj":r["pvals_adj"][g]}))
DE=pd.concat(out,ignore_index=True)
DE.to_csv(os.path.join(C.DATA,"LS_DE_wilcoxon.csv.gz"),index=False,compression="gzip")
up=DE[(DE.padj<0.001)&(DE.log2FC>1)].sort_values(["LS","log2FC"],ascending=[True,False])
up.to_csv(os.path.join(C.DATA,"LS_upgenes_wilcoxon.csv"),index=False)
print("wrote LS_upgenes_wilcoxon.csv; up-genes/state:\n",up.LS.value_counts().sort_index().to_string(),flush=True)
