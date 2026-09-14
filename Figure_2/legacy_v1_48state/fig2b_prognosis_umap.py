#!/usr/bin/env python
"""Figure 2B (v2): leukemic states coloured by prognosis (favorable/poor/ns) on the
full-atlas UMAP, same coordinates as Figure 1/2A. minou teal/red + grey (ns)."""
import os, numpy as np, pandas as pd, h5py, anndata as ad, scanpy as sc, matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
SC="/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
OUT="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_new_panels"
H5="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
f=h5py.File(H5,"r")
def col(n):
    nd=f["obs"][n]; cats=[c.decode() if isinstance(c,bytes) else c for c in nd["categories"][:]]
    return np.array([cats[i] if i>=0 else "NA" for i in nd["codes"][:]],dtype=object)
umap=f["obsm/X_umap"][:]; sub=col("AML Sub-Clusters"); mask=sub!="NA"; N=len(sub)
ci=col("clusters_for_interaction")
Xpca=f["obsm"]["X_pca"][:][:, :50].astype(np.float32); f.close()
# 48-state labels per leukemic cell (core + nearest-centroid)
Xl=Xpca[mask]; ref=pd.read_csv(f"{SC}/jaccard_hi/reference_assignments.csv")["res4.50"].astype(int).values
jc=pd.read_csv(f"{SC}/jaccard_hi/jaccard_percluster.csv"); jc=jc[jc.res==4.5]
stab=jc[jc.mean_jaccard>=0.35].sort_values("n_cells",ascending=False); raws=stab.cluster.astype(int).tolist()
raw2ls={rw:i for i,rw in enumerate(raws)}; is_core=np.isin(ref,raws)
cent=np.vstack([Xl[ref==rw].mean(0) for rw in raws])
lsidx=np.empty(len(ref),dtype=int); lsidx[is_core]=[raw2ls[r] for r in ref[is_core]]
lsidx[~is_core]=cdist(Xl[~is_core],cent.astype(np.float32)).argmin(1)
# state -> prognosis
pr=pd.read_csv(f"{SC}/target_ls/LS48_prognosis.csv")
prog_dir={int(r.LS.split("_")[1])-1:(r.dir if r.fdr<0.10 else "ns") for _,r in pr.iterrows()}
cellprog=np.array([prog_dir.get(i,"ns") for i in lsidx],dtype=object)
state=np.array([None]*N,dtype=object); state[np.where(mask)[0]]=cellprog
state[ci=="0_HSPC"]="normal HSPC"   # normal HSPC reference
adata=ad.AnnData(X=np.zeros((N,1),dtype=np.float32),
                 obs=pd.DataFrame({"Prognosis":pd.Categorical(state,categories=["favorable","poor","normal HSPC","ns"])}),
                 var=pd.DataFrame(index=["_dummy"]))
adata.obsm["X_umap"]=umap
adata.uns["Prognosis_colors"]=["#00798c","#d1495b","#E69F00","#d9d9d9"]  # favorable, poor, normal HSPC, ns
sc.settings.figdir=OUT; sc.set_figure_params(dpi=150,dpi_save=300,frameon=False)
sc.pl.umap(adata,color="Prognosis",title="Prognosis",frameon=False,show=False,
           groups=["favorable","poor","normal HSPC"],na_color="#e8e8e8")
plt.savefig(os.path.join(OUT,"Figure_2B_prognosis_umap.pdf"),bbox_inches="tight")
plt.savefig(os.path.join(OUT,"Figure_2B_prognosis_umap.png"),bbox_inches="tight",dpi=200); plt.close()
print("wrote Figure_2B_prognosis_umap; favorable %d poor %d ns %d cells"%(
    (cellprog=="favorable").sum(),(cellprog=="poor").sum(),(cellprog=="ns").sum()))
