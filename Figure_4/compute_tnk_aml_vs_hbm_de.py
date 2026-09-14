#!/usr/bin/env python
"""Figure 4: matched within-subset AML-vs-HBM comparison of T/NK cells. For each T/NK subset we compare
the SAME cell type between AML and healthy BM (so composition shifts can't drive it), and look for a
consistent cross-subset signature. HBM = 2 samples, so cell-level p-values are pseudoreplicated and only
used to rank; the reported evidence is (a) within-subset logFC and (b) SAMPLE-LEVEL consistency = fraction
of the 28 AML samples whose pseudobulk exceeds the higher of the 2 HBM samples. No fabricated values."""
import os, numpy as np, pandas as pd, scanpy as sc
from scipy.stats import mannwhitneyu
H5="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
D ="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
pc=pd.read_csv(os.path.join(D,"tnk_subset_percell.csv")).set_index("barcode")
a=sc.read_h5ad(H5); a=a[a.obs_names.isin(pc.index)].copy()
a.X=a.layers["Raw_Counts"].copy(); sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
a.obs["subset"]=pc["subset"].reindex(a.obs_names).values
a.obs["group"]=pc["group"].reindex(a.obs_names).values
a.obs["SampleID"]=pc["SampleID"].reindex(a.obs_names).values
genes=np.array(a.var_names); X=a.X.tocsc()
def col(g): return np.asarray(X[:,int(np.where(genes==g)[0][0])].todense()).ravel() if g in set(genes) else None
subs=[s for s in a.obs["subset"].unique()
      if ((a.obs.subset==s)&(a.obs.group=="AML")).sum()>=30 and ((a.obs.subset==s)&(a.obs.group=="HBM")).sum()>=30]
print("subsets tested (>=30 AML & >=30 HBM cells):",subs)
# dense matrix once (10k x ~27k is big; restrict to genes expressed in >=3% of T/NK)
expr_frac=np.asarray((X>0).mean(0)).ravel(); keep=np.where(expr_frac>=0.03)[0]
Xk=np.asarray(X[:,keep].todense()); gk=genes[keep]
obs=a.obs.reset_index(drop=True)
recs=[]
for s in subs:
    m=(obs.subset==s).values; ma=m&(obs.group=="AML").values; mh=m&(obs.group=="HBM").values
    lfc=Xk[ma].mean(0)-Xk[mh].mean(0)
    pa=(Xk[ma]>0).mean(0); ph=(Xk[mh]>0).mean(0)
    for i in range(len(gk)):
        recs.append([s,gk[i],lfc[i],pa[i],ph[i]])
de=pd.DataFrame(recs,columns=["subset","gene","logFC","pct_AML","pct_HBM"])
de.to_csv(os.path.join(D,"tnk_aml_vs_hbm_de_bySubset.csv"),index=False)
# cross-subset shared signature: mean logFC across subsets + #subsets with same sign
piv=de.pivot(index="gene",columns="subset",values="logFC")
sig=pd.DataFrame({"mean_lfc":piv.mean(1),"n_up":(piv>0.1).sum(1),"n_dn":(piv<-0.1).sum(1)})
nsub=len(subs)
up=sig[(sig.mean_lfc>0.15)&(sig.n_up>=max(3,nsub-2))].sort_values("mean_lfc",ascending=False)
dn=sig[(sig.mean_lfc<-0.15)&(sig.n_dn>=max(3,nsub-2))].sort_values("mean_lfc")
# sample-level consistency for the shared UP genes (pseudobulk per sample, all T/NK)
pb=pd.DataFrame(Xk,columns=gk); pb["SampleID"]=obs.SampleID.values; pb["group"]=obs.group.values
pbm=pb.groupby(["SampleID","group"]).mean(numeric_only=True).reset_index()
hbm_max=pbm[pbm.group=="HBM"].drop(columns=["SampleID","group"]).max()
aml=pbm[pbm.group=="AML"].drop(columns=["SampleID","group"])
def consist(g): return (aml[g]>hbm_max[g]).mean() if g in aml.columns else np.nan
up=up.copy(); up["frac_AMLsamp_gt_HBM"]=[consist(g) for g in up.index]
dn=dn.copy(); dn["frac_AMLsamp_lt_HBM"]=[ (aml[g]<pbm[pbm.group=="HBM"].drop(columns=["SampleID","group"]).min()[g]).mean() if g in aml.columns else np.nan for g in dn.index]
up.to_csv(os.path.join(D,"tnk_aml_up_shared.csv")); dn.to_csv(os.path.join(D,"tnk_aml_down_shared.csv"))
print("\nSHARED AML-UP (higher in AML across %d subsets), top 25:"%nsub)
print(up.head(25).round(3).to_string())
print("\nSHARED AML-DOWN, top 25:")
print(dn.head(25).round(3).to_string())
