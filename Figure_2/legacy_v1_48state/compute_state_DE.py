#!/usr/bin/env python
"""Per-cell Wilcoxon DE for the NEW LS states (state vs rest), matching the
original Figure 2F method. Core cells only (cells in the stable res-4.5 cluster).
Outputs full DE table + up-gene lists (padj<0.001 & log2FC>1)."""
import scanpy as sc, numpy as np, pandas as pd, h5py, scipy.sparse as sp
SC="/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
H5="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"

# --- state labels per leukemic cell (core = cell's res4.5 raw is a stable cluster) ---
f=h5py.File(H5,"r")
nd=f["obs"]["AML Sub-Clusters"]
cats=[c.decode() if isinstance(c,bytes) else c for c in nd["categories"][:]]
sub=np.array([cats[i] if i>=0 else "NA" for i in nd["codes"][:]],dtype=object)
f.close()
mask=sub!="NA"; N=len(sub); pos=np.where(mask)[0]
ref=pd.read_csv(f"{SC}/jaccard_hi/reference_assignments.csv")["res4.50"].astype(int).values
jac=pd.read_csv(f"{SC}/jaccard_hi/jaccard_percluster.csv"); jac=jac[jac.res==4.5]
stab=jac[jac.mean_jaccard>=0.35].sort_values("n_cells",ascending=False)
raws=stab.cluster.astype(int).tolist(); ls=[f"LS_{i+1}" for i in range(len(raws))]
raw2ls={rw:i for i,rw in enumerate(raws)}
is_core=np.isin(ref,raws)
keep=pos[is_core]
lslab=[ls[raw2ls[r]] for r in ref[is_core]]
print(f"leukemic={mask.sum()}, core-labelled={len(keep)}, states={len(ls)}",flush=True)

# --- load only the core cells (backed subset), ensure log-normalized X ---
adata=sc.read_h5ad(H5,backed="r")
ad=adata[keep].to_memory()
ad.obs["LS"]=pd.Categorical(lslab,categories=ls)
mn=float(ad.X.min()); mx=float(ad.X.max())
print(f"X min={mn:.3f} max={mx:.3f}",flush=True)
if mn<-0.01:                                   # scaled -> rebuild from counts
    lyr="counts" if "counts" in ad.layers else ("Raw_Counts" if "Raw_Counts" in ad.layers else None)
    ad.X=ad.layers[lyr].copy(); sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad); print("rebuilt lognorm from",lyr,flush=True)
elif mx>50:                                    # raw counts -> normalize+log
    sc.pp.normalize_total(ad,target_sum=1e4); sc.pp.log1p(ad); print("normalized+log1p from counts-like X",flush=True)
else:
    print("X already log-normalized",flush=True)

sc.tl.rank_genes_groups(ad,"LS",method="wilcoxon")
r=ad.uns["rank_genes_groups"]; groups=r["names"].dtype.names
out=[]
for g in groups:
    out.append(pd.DataFrame({"LS":g,"gene":r["names"][g],
                             "log2FC":r["logfoldchanges"][g],"padj":r["pvals_adj"][g]}))
DE=pd.concat(out,ignore_index=True)
DE.to_csv(f"{SC}/ls48/LS48_DE_wilcoxon.csv.gz",index=False,compression="gzip")
up=DE[(DE.padj<0.001)&(DE.log2FC>1)].copy().sort_values(["LS","log2FC"],ascending=[True,False])
up.to_csv(f"{SC}/ls48/LS48_upgenes_wilcoxon.csv",index=False)
print("up-genes per state (padj<0.001 & log2FC>1):",flush=True)
print(up.LS.value_counts().sort_index().to_string(),flush=True)
print("DONE",flush=True)
