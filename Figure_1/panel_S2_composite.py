#!/usr/bin/env python
"""S2 (AML atlas covariates), rebuilt for the corrected scheme. (A) atlas UMAP by sample;
(B) covariate UMAPs: age, outcome, RELAPSE (authoritative clinical annotation), remission;
(C) pseudotime (dpt, palantir, entropy); (D) lineage-marker dot plot per cell type. Relapse
uses the reconciled clinical set (memory relapse-annotation-source; AML3121 corrected), not the
unreliable h5ad field. Data = h5ad obs/obsm + 'counts' layer. No fabricated values."""
import os, numpy as np, pandas as pd, h5py, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.sparse import csr_matrix
H5="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
OUT="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/05_Supplementary_Figures"
RELAPSED={"AML948","AML3082","AML3371","AML3492","AML4068","AML4239","AML4304","AML3121"}
UNKNOWN={"AML4271","AML4363"}
MARK=["PTPRC","CD34","AVP","CD38","MPO","ELANE","AZU1","LYZ","CD14","FCN1","C1QA","CD1C","FCER1A",
      "IRF7","LILRA4","GATA1","HBB","GYPA","PF4","ITGA2B","CD19","MS4A1","CD79A","CD3D","IL7R","NKG7","GNLY","MKI67","PCNA","CD96"]
CTORD=["AML","AML-MKI67","AML-PCNA","AML-CD1C","AML-CD14","AML-B","AML-Ery","AML-NK","AML-CTL","AML-CD4T","AML-Naïve_CD8T",
       "0_HSPC","Myeloid_Pro","CD14_Monocyte","CD16_Monocyte","Macrophage","mDC","pDC","Erythrocytes",
       "CD34+ProB","ProB","PreB","CD20+B","PlasmaB","Naïve_CD4T","Activated_CD4T","Naïve_CD8T","CTL","NK"]

f=h5py.File(H5,"r")
def dec(a): return np.array([x.decode() if isinstance(x,bytes) else x for x in a])
def cat(name):
    g=f["obs"][name]; return pd.Series(dec(g["categories"][:])[g["codes"][:]])
obs=dec(f["obs"]["_index"][:]); var=dec(f["var"]["_index"][:])
um=f["obsm"]["X_umap"][:]
df=pd.DataFrame({"sample":cat("SampleID").values,"age":cat("Age").values,"prognosis":cat("Prognosis").values,
                "remission":cat("Remission").values,"celltype":cat("Cell_Type").values,
                "dpt":f["obs"]["dpt_pseudotime"][:],"pal":f["obs"]["palantir_pseudotime"][:],"ent":f["obs"]["palantir_entropy"][:]}, index=obs)
df["U1"]=um[:,0]; df["U2"]=um[:,1]
def relc(s): return "healthy BM" if s.startswith("0_HealthyBM") else ("relapsed" if s in RELAPSED else ("unknown" if s in UNKNOWN else "not relapsed"))
df["relapse"]=[relc(s) for s in df["sample"]]
# markers expression (CSR counts, log-norm)
gpos={g:i for i,g in enumerate(var)}; gi=[gpos[g] for g in MARK if g in gpos]; present=[g for g in MARK if g in gpos]
grp=f["layers"]["counts"]; shape=tuple(grp.attrs["shape"])
X=csr_matrix((grp["data"][:],grp["indices"][:],grp["indptr"][:]),shape=shape)
f.close()
E=pd.DataFrame(np.asarray(X[:,gi].todense()),columns=present,index=obs); E["celltype"]=df["celltype"].values

fig=plt.figure(figsize=(15,19)); gs=GridSpec(4,4,figure=fig,height_ratios=[1.25,1,1,1.5],hspace=0.35,wspace=0.28)
def umap_cat(ax,key,title,cmap=None,order=None,leg=True,s=1.4):
    vals=df[key].astype(str); cats=order or sorted(vals.unique())
    cols=cmap or {c:plt.cm.tab20(i%20) for i,c in enumerate(cats)}
    for c in cats:
        m=vals==c; ax.scatter(df.U1[m],df.U2[m],s=s,color=cols[c],linewidths=0,rasterized=True,label=c)
    ax.set_xticks([]);ax.set_yticks([]);ax.set_title(title,fontsize=10)
    if leg: ax.legend(loc="upper right",fontsize=5.5,markerscale=4,frameon=False,handletextpad=0.1)
def umap_cont(ax,key,title):
    sc=ax.scatter(df.U1,df.U2,s=1.4,c=df[key],cmap="viridis",linewidths=0,rasterized=True)
    ax.set_xticks([]);ax.set_yticks([]);ax.set_title(title,fontsize=10); fig.colorbar(sc,ax=ax,shrink=0.6)
# A sample (span 2 cols)
axA=fig.add_subplot(gs[0,0:2]); umap_cat(axA,"sample","A  Sample",leg=True,s=1.2)
# B covariates
relcol={"not relapsed":"#c7ccd1","relapsed":"#d1495b","unknown":"#8d99ae","healthy BM":"#2a9d8f"}
progcol={"Alive":"#4a6fe3","Deceased":"#d1495b","0_HealthyBM":"#2a9d8f"}
remcol={"True":"#4393c3","False":"#d6604d","0_HealthyBM":"#2a9d8f"}
agecol={"0-1":"#fee5d9","1-10":"#fcae91","10":"#fb6a4a",">10":"#de2d26",">20":"#a50f15","0_HealthyBM":"#2a9d8f"}
umap_cat(fig.add_subplot(gs[0,2]),"age","B  Age",cmap=agecol,order=["0-1","1-10",">10",">20","0_HealthyBM"],s=1.2)
umap_cat(fig.add_subplot(gs[0,3]),"prognosis","Outcome",cmap=progcol,order=["Alive","Deceased","0_HealthyBM"],s=1.2)
umap_cat(fig.add_subplot(gs[1,0]),"relapse","Relapse (clinical, corrected)",cmap=relcol,order=["not relapsed","relapsed","unknown","healthy BM"],s=1.2)
umap_cat(fig.add_subplot(gs[1,1]),"remission","Remission",cmap=remcol,order=["True","False","0_HealthyBM"],s=1.2)
# C pseudotime
umap_cont(fig.add_subplot(gs[1,2]),"dpt","C  dpt pseudotime")
umap_cont(fig.add_subplot(gs[1,3]),"pal","palantir pseudotime")
umap_cont(fig.add_subplot(gs[2,0]),"ent","palantir entropy")
# D dot plot (markers x cell type)
axD=fig.add_subplot(gs[3,0:4])
cts=[c for c in CTORD if c in set(E["celltype"])]
mean=E.groupby("celltype")[present].mean().reindex(cts)
pct=E.groupby("celltype")[present].apply(lambda d:(d>0).mean()).reindex(cts)
z=(mean-mean.mean(0))/(mean.std(0)+1e-9)
for yi,ct in enumerate(cts):
    for xi,g in enumerate(present):
        axD.scatter(xi,yi,s=pct.loc[ct,g]*90+1,c=[z.loc[ct,g]],cmap="Reds",vmin=-1,vmax=2.5,edgecolors="none")
axD.set_xticks(range(len(present))); axD.set_xticklabels(present,rotation=90,fontsize=7)
axD.set_yticks(range(len(cts))); axD.set_yticklabels(cts,fontsize=7); axD.invert_yaxis()
axD.set_title("D  Lineage-marker expression per cell type (size = % expressing, colour = z mean)",fontsize=10)
axD.set_xlim(-1,len(present)); axD.set_ylim(len(cts),-1)
for e in ("png","pdf"): fig.savefig(os.path.join(OUT,f"Figure_S2__AML_atlas_covariates.{e}"),dpi=140,bbox_inches="tight")
print("wrote rebuilt S2 (relapsed cells=%d)"%(df['relapse']=='relapsed').sum())
