#!/usr/bin/env python
"""Figure 5 panel D redo on the corrected scheme: candidate-target coverage (% positive + mean log-norm)
across the 8 normal hematopoietic compartments (toxicity reference) and the leukemic compartment
(efficacy), replacing the retired PAC/PPAC/FPAC populations. Gene symbols mapped: CD123->IL3RA,
CLL-1->CLEC12A. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
sys.path.insert(0,"/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"); import config as C
OUT="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_5"
GSET=["CD96","CD9","SUCNR1","AMN","IL1RAP","LMBR1L","NRXN2","FURIN","UMODL1","CD7","CD33","CD123","CLL-1","CD38","CD70","FLT3","TNFRSF4","ABCA7","ITGAX"]
SYM={"CD123":"IL3RA","CLL-1":"CLEC12A"}   # figure label -> gene symbol
NORM={"0_HSPC":"HSPC","Myeloid_Pro":"Myeloid_Pro","CD14_Monocyte":"Monocyte","CD16_Monocyte":"Monocyte",
      "Macrophage":"Macrophage","mDC":"DC","pDC":"DC","Erythrocytes":"Erythroid",
      "Naïve_CD4T":"T/NK","Naïve_CD8T":"T/NK","CTL":"T/NK","NK":"T/NK","Activated_CD4T":"T/NK",
      "CD20+B":"B","ProB":"B","PreB":"B","PlasmaB":"B","CD34+ProB":"B"}
f=h5py.File(C.H5,"r"); dec=lambda a:np.array([x.decode() if isinstance(x,bytes) else x for x in a])
obs=dec(f["obs"]["_index"][:]); var=dec(f["var"]["_index"][:])
ct=pd.Series(dec(f["obs"]["Cell_Type"]["categories"][:])[f["obs"]["Cell_Type"]["codes"][:]],index=obs)
st=pd.Series(dec(f["obs"]["SampleType"]["categories"][:])[f["obs"]["SampleType"]["codes"][:]],index=obs)
gpos={g:i for i,g in enumerate(var)}
present=[(lab, SYM.get(lab,lab)) for lab in GSET if SYM.get(lab,lab) in gpos]
gi=[gpos[sym] for _,sym in present]
grp=f["layers"]["counts"]; X=csr_matrix((grp["data"][:],grp["indices"][:],grp["indptr"][:]),shape=tuple(grp.attrs["shape"])); f.close()
E=pd.DataFrame(np.asarray(X[:,gi].todense()),columns=[lab for lab,_ in present],index=obs)
# leukemic = AML-sample cells in the leukemic Cell_Type set
LEUK={"AML","AML-MKI67","AML-PCNA","AML-CD1C"}
grpcol=pd.Series(index=obs,dtype=object)
grpcol[(ct.isin(LEUK)) & (st=="AML")]="leukemic"
for k,v in NORM.items(): grpcol[ct==k]=v
E["group"]=grpcol.values; E=E[E["group"].notna()]
order=["leukemic","HSPC","Myeloid_Pro","Monocyte","Macrophage","DC","Erythroid","B","T/NK"]
rows=[]
for lab,_ in present:
    for g in order:
        d=E.loc[E["group"]==g,lab]
        if len(d)==0: continue
        rows.append([lab,g,round(100*(d>0).mean(),2),round(d.mean(),4)])
out=pd.DataFrame(rows,columns=["gene","group","pct","mean_expr"]); out["label"]=out["gene"]
out.to_csv(os.path.join(OUT,"fig5D_ls.csv"),index=False)
print("wrote fig5D_ls.csv | groups:", order)
print("\nCD96 by group (pct positive / mean):")
print(out[out.gene=="CD96"][["group","pct","mean_expr"]].to_string(index=False))
