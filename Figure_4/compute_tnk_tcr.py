#!/usr/bin/env python
"""Figure 4 TCR analysis. From the scirpy-format T-cell h5ad (paired TRA/TRB), define clonotypes by
productive CDR3 amino-acid sequence kept PATIENT-PRIVATE (keyed on SampleID, since TCRs are patient
specific), compute clonal expansion and a clonality index (1 - normalized Shannon entropy) per T/NK
subset and by AML vs healthy BM, and test whether clonally-expanded / exhausted T cells co-express CD96
(the lead leukemic target, also an inhibitory receptor in the TIGIT/CD226/CD96 axis). No fabrication."""
import os, numpy as np, pandas as pd, scanpy as sc
H5="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_total_Tcell_rna_withTCR.h5ad"
D ="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
CKPT=["CD96","TIGIT","CD226","PDCD1","HAVCR2","LAG3","TOX"]
a=sc.read_h5ad(H5)
print("X max (expect log-norm ~<10):", round(float(a.X.max()),2))
sub=pd.read_csv(os.path.join(D,"tnk_subset_percell.csv")).set_index("barcode")
a=a[a.obs_names.isin(sub.index)].copy()
o=a.obs
o["subset"]=sub["subset"].reindex(a.obs_names).values
o["group"]=sub["group"].reindex(a.obs_names).values
# --- clonotype: productive TRB (required) + TRA if present, kept patient-private ---
def prod(col,pcol):
    s=o[col].astype(str); p=o[pcol].astype(str)
    s=s.where((p=="True") & s.notna() & (s!="None") & (s!="nan"), other="")
    return s
trb=prod("IR_VDJ_1_junction_aa","IR_VDJ_1_productive")
tra=prod("IR_VJ_1_junction_aa","IR_VJ_1_productive")
o["has_tcr"]=trb!=""
o["clonotype"]=np.where(o["has_tcr"], o["SampleID"].astype(str)+"|"+tra+"_"+trb, np.nan)
tc=o[o["has_tcr"]].copy()
tc["clone_size"]=tc.groupby("clonotype")["clonotype"].transform("size")
tc["expanded"]=np.where(tc["clone_size"]>=2,"expanded","unexpanded")
o.loc[tc.index,"clone_size"]=tc["clone_size"]; o.loc[tc.index,"expanded"]=tc["expanded"]
print("cells with productive TRB:",int(o["has_tcr"].sum()),"| unique clonotypes:",tc["clonotype"].nunique())

def clonality(sizes):
    # sizes = per-clonotype cell counts; clonality = 1 - Shannon/ln(#clones)
    n=len(sizes)
    if n<=1: return np.nan
    p=np.asarray(sizes,float); p=p/p.sum()
    H=-(p*np.log(p)).sum()
    return 1 - H/np.log(n)

# per-subset expansion + clonality
rows=[]
for s,g in tc.groupby("subset"):
    sizes=g.groupby("clonotype").size().values
    rows.append([s,len(g),(g["clone_size"]>=2).mean(),clonality(sizes)])
bysub=pd.DataFrame(rows,columns=["subset","n_tcr","frac_expanded","clonality"]).sort_values("clonality",ascending=False)
bysub.to_csv(os.path.join(D,"tcr_by_subset.csv"),index=False)
print("\n== expansion + clonality by subset ==\n",bysub.round(3).to_string(index=False))

# per-sample clonality, AML vs HBM
rows=[]
for (samp,grp),g in tc.groupby(["SampleID","group"]):
    if len(g)<20: continue
    rows.append([samp,grp,len(g),clonality(g.groupby("clonotype").size().values),(g["clone_size"]>=2).mean()])
bysamp=pd.DataFrame(rows,columns=["SampleID","group","n_tcr","clonality","frac_expanded"])
bysamp.to_csv(os.path.join(D,"tcr_clonality_by_sample.csv"),index=False)
print("\n== clonality AML vs HBM (per-sample medians) ==")
print(bysamp.groupby("group")[["clonality","frac_expanded","n_tcr"]].median().round(3).to_string())

# --- CD96 / checkpoint co-expression from the MAIN atlas counts layer (log-norm); the withTCR h5ad X is
#     scaled/z-scored (max ~86, negatives) so is unusable for expression magnitude ---
import h5py
from scipy.sparse import csr_matrix
MAIN="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
f=h5py.File(MAIN,"r"); dec=lambda x:[v.decode() if isinstance(v,bytes) else v for v in x]
mobs=np.array(dec(f["obs"]["_index"][:])); mvar=np.array(dec(f["var"]["_index"][:]))
present=[g for g in CKPT if g in set(mvar)]
gi=[int(np.where(mvar==g)[0][0]) for g in present]
gg=f["layers"]["counts"]; X=csr_matrix((gg["data"][:],gg["indices"][:],gg["indptr"][:]),shape=tuple(gg.attrs["shape"]))
f.close()
Em=pd.DataFrame(np.asarray(X[:,gi].todense()),columns=present,index=mobs)
print("main-atlas counts layer max (expect log-norm):", round(float(Em.values.max()),2))
E=Em.reindex(a.obs_names).copy()
E["subset"]=o["subset"].values; E["expanded"]=o["expanded"].values; E["group"]=o["group"].values
# by subset (mean + pct)
rows=[]
for s,g in E.groupby("subset"):
    for gene in present: rows.append([s,gene,g[gene].mean(),(g[gene]>0).mean()])
pd.DataFrame(rows,columns=["subset","gene","mean","pct"]).to_csv(os.path.join(D,"tcr_ckpt_by_subset.csv"),index=False)
# expanded vs unexpanded within cytotoxic CD8 (GZMK+GZMB CD8 T)
cd8=E[E.subset.isin(["GZMK CD8 T","GZMB CD8 T","Memory CD8 T"]) & E.expanded.notna()]
print("\n== checkpoint mean log-norm in cytotoxic/memory CD8: expanded vs unexpanded ==")
print(cd8.groupby("expanded")[present].mean().round(3).to_string())
rows=[]
for (s,ex),g in E[E.expanded.notna()].groupby(["subset","expanded"]):
    for gene in present: rows.append([s,ex,gene,g[gene].mean(),(g[gene]>0).mean()])
pd.DataFrame(rows,columns=["subset","expanded","gene","mean","pct"]).to_csv(os.path.join(D,"tcr_ckpt_by_expansion.csv"),index=False)
print("\nwrote tcr_by_subset.csv / tcr_clonality_by_sample.csv / tcr_ckpt_by_subset.csv / tcr_ckpt_by_expansion.csv")
