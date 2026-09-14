#!/usr/bin/env python
"""Figure 2A (v2), rendered exactly like Figure 1 (rebuild_figure1_umaps.py):
sc.pl.umap on the FULL atlas from obsm/X_umap, frameon=False. Two panels:
CD34 expression, and the 30 leukemic states (non-leukemic cells greyed as NA).
Same coordinates as Figure 1 -> no need to touch Figure 1."""
import os, numpy as np, pandas as pd, h5py, anndata as ad, scanpy as sc, matplotlib, matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.spatial.distance import cdist
from matplotlib.colors import LinearSegmentedColormap
HEATMAP0=LinearSegmentedColormap.from_list("heatmap0",["#001219","#005F73","#0A9396","#94D2BD","#E9D8A6","#EE9B00","#CA6702","#AE2012","#9B2226"])
SC="/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
OUT="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_new_panels"
H5="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
NSTATE=48
pal30=[matplotlib.colors.to_hex(c) for c in (list(plt.cm.tab20.colors)+list(plt.cm.tab20b.colors)+list(plt.cm.tab20c.colors))][:NSTATE]

f=h5py.File(H5,"r")
def col(n):
    nd=f["obs"][n]
    if isinstance(nd,h5py.Group):
        cats=[c.decode() if isinstance(c,bytes) else c for c in nd["categories"][:]]; return np.array([cats[i] if i>=0 else "NA" for i in nd["codes"][:]],dtype=object)
    a=nd[:]; return np.array([x.decode() for x in a],dtype=object) if a.dtype.kind=="S" else a
umap=f["obsm/X_umap"][:]; sub=col("AML Sub-Clusters"); mask=sub!="NA"; N=len(sub)
Xpca=f["obsm"]["X_pca"][:][:, :50].astype(np.float32)
genes=np.array([g.decode() if isinstance(g,bytes) else g for g in f["var"]["_index"][:]])
cd34=int(np.where(genes=="CD34")[0][0])
RC=f["layers"]["Raw_Counts"]; dd,ii,ip=RC["data"],RC["indices"],RC["indptr"][:]; G=len(genes)
cvec=np.zeros(N,dtype=np.float32)
for a in range(0,N,8000):
    b=min(a+8000,N); s,e=int(ip[a]),int(ip[b])
    X=csr_matrix((dd[s:e].astype(np.float32),ii[s:e],(ip[a:b+1]-ip[a]).astype(np.int64)),shape=(b-a,G))
    rs=np.asarray(X.sum(1)).ravel(); rs[rs==0]=1
    cvec[a:b]=np.asarray(X[:,cd34].todense()).ravel()*(1e6/rs)
f.close()

# 30-state labels for leukemic cells (core + nearest-centroid), NA elsewhere
Xl=Xpca[mask]; ref=pd.read_csv(f"{SC}/jaccard_hi/reference_assignments.csv")["res4.50"].astype(int).values
pc=pd.read_csv(f"{SC}/jaccard_hi/jaccard_percluster.csv"); pc=pc[pc.res==4.5]
stab=pc[pc.mean_jaccard>=0.35].sort_values("n_cells",ascending=False); raws=stab.cluster.astype(int).tolist()
raw2ls={rw:i for i,rw in enumerate(raws)}; is_core=np.isin(ref,raws)
cent=np.vstack([Xl[ref==rw].mean(0) for rw in raws])
lsidx=np.empty(len(ref),dtype=int); lsidx[is_core]=[raw2ls[r] for r in ref[is_core]]
lsidx[~is_core]=cdist(Xl[~is_core],cent.astype(np.float32)).argmin(1)
state=np.array([None]*N,dtype=object); state[np.where(mask)[0]]=[f"LS{i+1}" for i in lsidx]

adata=ad.AnnData(X=np.zeros((N,1),dtype=np.float32),
                 obs=pd.DataFrame({"Leukemic states":pd.Categorical(state,categories=[f"LS{i+1}" for i in range(NSTATE)]),
                                   "CD34":np.log1p(cvec)}),
                 var=pd.DataFrame(index=["_dummy"]))
adata.obsm["X_umap"]=umap; adata.uns["Leukemic states_colors"]=pal30
sc.settings.figdir=OUT; sc.set_figure_params(dpi=150,dpi_save=300,frameon=False)
sc.pl.umap(adata,color="CD34",title="CD34",frameon=False,show=False,cmap=HEATMAP0)
plt.savefig(os.path.join(OUT,"Figure_2A_CD34.pdf"),bbox_inches="tight"); plt.savefig(os.path.join(OUT,"Figure_2A_CD34.png"),bbox_inches="tight",dpi=200); plt.close()
sc.pl.umap(adata,color="Leukemic states",title="48 leukemic states",frameon=False,show=False,legend_loc="on data",legend_fontsize=6,legend_fontoutline=2)
plt.savefig(os.path.join(OUT,"Figure_2A_states.pdf"),bbox_inches="tight"); plt.savefig(os.path.join(OUT,"Figure_2A_states.png"),bbox_inches="tight",dpi=200); plt.close()
sc.pl.umap(adata,color=["CD34","Leukemic states"],frameon=False,show=False,wspace=0.3,cmap=HEATMAP0)
plt.savefig(os.path.join(OUT,"Figure_2A_combined.png"),bbox_inches="tight",dpi=200); plt.close()
print("wrote Figure_2A_* (scanpy, full-atlas, Fig1-consistent) to",OUT)
