#!/usr/bin/env python
"""Figure 4 extended TCR analyses:
 (1) Shannon-based clonality per subset, DOWNSAMPLED to equal n (bootstrap) for a fair comparison,
     plus Inverse-Simpson as a cross-check.
 (2) Clonal SHARING between subsets (patient-private clonotypes) - TCR as a lineage barcode
     (esp. GZMK-CD8 <-> GZMB-CD8 differentiation).
 (3) Are big clones the exhausted ones - exhaustion score (PDCD1/HAVCR2/LAG3/TOX, main-atlas log-norm)
     vs clonal expansion.
 (4) Clonality by cytogenetic subtype (per sample).
 (5) TRAV / TRBV V-gene usage, AML vs HBM (the original panel D).
No fabricated values."""
import os, numpy as np, pandas as pd, scanpy as sc, h5py
from scipy.sparse import csr_matrix
H5="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_total_Tcell_rna_withTCR.h5ad"
MAIN="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
D ="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
rng=np.random.RandomState(0)
a=sc.read_h5ad(H5)
sub=pd.read_csv(os.path.join(D,"tnk_subset_percell.csv")).set_index("barcode")
a=a[a.obs_names.isin(sub.index)].copy(); o=a.obs
o["subset"]=sub["subset"].reindex(a.obs_names).values
o["group"]=sub["group"].reindex(a.obs_names).values
def prod(col,pcol):
    s=o[col].astype(str); p=o[pcol].astype(str)
    return s.where((p=="True")&s.notna()&(s!="None")&(s!="nan"), other="")
trb=prod("IR_VDJ_1_junction_aa","IR_VDJ_1_productive"); tra=prod("IR_VJ_1_junction_aa","IR_VJ_1_productive")
o["has_tcr"]=trb!=""
o["clonotype"]=np.where(o["has_tcr"], o["SampleID"].astype(str)+"|"+tra+"_"+trb, np.nan)
o["trav"]=prod("IR_VJ_1_v_call","IR_VJ_1_productive"); o["trbv"]=prod("IR_VDJ_1_v_call","IR_VDJ_1_productive")
tc=o[o["has_tcr"]].copy()
tc["clone_size"]=tc.groupby("clonotype")["clonotype"].transform("size")

def clonality(sizes):
    n=len(sizes)
    if n<=1: return np.nan
    p=np.asarray(sizes,float); p=p/p.sum()
    return 1 - (-(p*np.log(p)).sum())/np.log(n)
def invsimpson(sizes):
    p=np.asarray(sizes,float); p=p/p.sum(); return 1/np.sum(p**2)

# (1) downsampled clonality per subset (n=100, bootstrap 200)
NMIN, NBOOT = 100, 200
counts=tc["subset"].value_counts()
elig=[s for s in counts.index if counts[s]>=NMIN]
rows=[]
for s in elig:
    idx=tc.index[tc["subset"]==s].to_numpy()
    cl=[]; iv=[]
    for _ in range(NBOOT):
        samp=rng.choice(idx, NMIN, replace=False)
        sizes=pd.Series(tc.loc[samp,"clonotype"]).value_counts().values
        cl.append(clonality(sizes)); iv.append(invsimpson(sizes))
    rows.append([s,counts[s],np.mean(cl),np.std(cl),np.mean(iv)])
pd.DataFrame(rows,columns=["subset","n_tcr","clonality_ds","clonality_sd","invsimpson_ds"]).sort_values("clonality_ds",ascending=False).to_csv(os.path.join(D,"tcr_clonality_downsampled.csv"),index=False)
print("(1) downsampled clonality written for",len(elig),"subsets (n>=%d)"%NMIN)

# (2) clonal sharing between subsets (Jaccard of clonotype sets)
sets={s:set(tc.loc[tc["subset"]==s,"clonotype"]) for s in counts.index if counts[s]>=20}
subs=list(sets); rows=[]
for i in subs:
    for j in subs:
        inter=len(sets[i]&sets[j]); uni=len(sets[i]|sets[j])
        rows.append([i,j,inter,inter/uni if uni else 0, inter/len(sets[i]) if len(sets[i]) else 0])
pd.DataFrame(rows,columns=["subset_a","subset_b","n_shared","jaccard","frac_a_in_b"]).to_csv(os.path.join(D,"tcr_sharing.csv"),index=False)
gk=len(sets.get("GZMK CD8 T",set())&sets.get("GZMB CD8 T",set()))
print("(2) sharing written | GZMK-CD8 & GZMB-CD8 share",gk,"clonotypes")

# (3) exhaustion score (main-atlas log-norm) vs clonal expansion
EXH=["PDCD1","HAVCR2","LAG3","TOX"]
f=h5py.File(MAIN,"r"); dec=lambda x:[v.decode() if isinstance(v,bytes) else v for v in x]
mobs=np.array(dec(f["obs"]["_index"][:])); mvar=np.array(dec(f["var"]["_index"][:]))
gi=[int(np.where(mvar==g)[0][0]) for g in EXH if g in set(mvar)]
gg=f["layers"]["counts"]; X=csr_matrix((gg["data"][:],gg["indices"][:],gg["indptr"][:]),shape=tuple(gg.attrs["shape"])); f.close()
exh=pd.Series(np.asarray(X[:,gi].todense()).mean(1).ravel(),index=mobs)
tc["exhaustion"]=exh.reindex(tc.index).values
tc["expanded"]=np.where(tc["clone_size"]>=2,"expanded","unexpanded")
ex=tc.groupby(["subset","expanded"])["exhaustion"].mean().unstack()
ex.to_csv(os.path.join(D,"tcr_exhaustion_by_expansion.csv"))
tc[["subset","group","clone_size","exhaustion","expanded"]].to_csv(os.path.join(D,"tcr_percell_exh.csv"))
print("(3) exhaustion vs expansion (cytotoxic CD8):")
print(ex.reindex(["GZMK CD8 T","GZMB CD8 T","Memory CD8 T"]).round(3).to_string())

# (4) clonality by cytogenetic subtype (per AML sample)
o2=tc.copy(); o2["Cytogenetic"]=a.obs["Cytogenetic"].reindex(tc.index).astype(str).values
rows=[]
for (samp,cyto),g in o2[o2["group"]=="AML"].groupby(["SampleID","Cytogenetic"]):
    if len(g)>=20: rows.append([samp,cyto,len(g),clonality(g.groupby("clonotype").size().values)])
pd.DataFrame(rows,columns=["SampleID","Cytogenetic","n_tcr","clonality"]).to_csv(os.path.join(D,"tcr_clonality_by_cyto.csv"),index=False)
print("(4) clonality by cytogenetic written")

# (5) TRBV / TRAV usage AML vs HBM (fraction of TCR+ cells per V gene)
for chain,coln,fn in [("TRBV","trbv","tcr_trbv_usage.csv"),("TRAV","trav","tcr_trav_usage.csv")]:
    d=tc[tc[coln]!=""].copy()
    u=d.groupby(["group",coln]).size().rename("n").reset_index()
    u["frac"]=u["n"]/u.groupby("group")["n"].transform("sum")
    u.rename(columns={coln:"Vgene"}).to_csv(os.path.join(D,fn),index=False)
print("(5) V-gene usage written")
