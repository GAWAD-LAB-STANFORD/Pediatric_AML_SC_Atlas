#!/usr/bin/env python
"""Figure 3G data: per-group % expressing + mean log-CPM for four program modules --
HSC quiescence (HLF/AVP/CRHBP), leukemic stemness (HOXA9/MEIS1), proliferation
(MKI67/TOP2A/CDK1), myeloid differentiation (MPO/LYZ/CD14) -- across normal HSPC, the
leukemic HSC/MPP states (LS5/LS38), LS10, and other leukemic cells. Reads raw counts
straight from the h5ad. No fabricated values."""
import numpy as np, pandas as pd, h5py
H5="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
D3="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3/"
OUT=D3+"LS_program_modules.csv"
MODULES={"HSC quiescence":["HLF","AVP","CRHBP"], "leukemic stemness":["HOXA9","MEIS1"],
         "proliferation":["MKI67","TOP2A","CDK1"], "differentiation":["MPO","LYZ","CD14"]}
MK=[g for gs in MODULES.values() for g in gs]
ra=pd.read_csv(D3+"reference_assignments.csv"); ra["barcode"]=ra["barcode"].astype(str).str.replace(r"^b'|'$","",regex=True)
jp=pd.read_csv(D3+"jaccard_percluster.csv"); jp=jp[jp.res==4.5]
raws=jp[jp.mean_jaccard>=0.35].sort_values("n_cells",ascending=False).cluster.astype(int).tolist()
raw2ls={rw:f"LS_{i+1}" for i,rw in enumerate(raws)}; ra["LS"]=ra["res4.50"].map(raw2ls)
bc2ls=dict(zip(ra.barcode,ra.LS))
f=h5py.File(H5,"r")
names=np.array([x.decode() if isinstance(x,bytes) else x for x in f["obs"]["_index"][:]])
nd=f["obs"]["Cell Type"]; cc=[x.decode() if isinstance(x,bytes) else x for x in nd["categories"][:]]
ctype=np.array([cc[i] if i>=0 else "NA" for i in nd["codes"][:]],dtype=object)
genes=np.array([g.decode() if isinstance(g,bytes) else g for g in f["var"]["_index"][:]])
gi={g:int(np.where(genes==g)[0][0]) for g in MK if g in set(genes)}
RC=f["layers"]["Raw_Counts"]; dd,ii,ip=RC["data"],RC["indices"],RC["indptr"][:]
def groupof(b,ct):
    if ct=="0_HSPC": return "normal HSPC"
    ls=bc2ls.get(b)
    if ls in ("LS_5","LS_38"): return "LS5/LS38 (leuk. HSC/MPP)"
    if ls=="LS_10": return "LS10 (proliferative)"
    if isinstance(ls,str): return "other leukemic"
    return None
grp=np.array([groupof(b,ct) for b,ct in zip(names,ctype)],dtype=object)
sel=np.where(grp!=None)[0]
E={g:np.zeros(len(sel)) for g in gi}
for k,c in enumerate(sel):
    s,e=int(ip[c]),int(ip[c+1]); ri=ii[s:e]; rd=dd[s:e]; t=rd.sum()
    for g,gg in gi.items():
        v=rd[ri==gg].sum() if (ri==gg).any() else 0; E[g][k]=np.log1p(v/t*1e6) if t>0 else 0
f.close()
d=pd.DataFrame(E); d["grp"]=grp[sel]
gene2mod={g:m for m,gs in MODULES.items() for g in gs}
rows=[]
for g in gi:
    for gr,x in d.groupby("grp"):
        rows.append([g,gene2mod[g],gr,round((x[g]>0).mean()*100,1),round(x[g].mean(),3),len(x)])
out=pd.DataFrame(rows,columns=["gene","module","group","pct","mean_expr","n"])
out.to_csv(OUT,index=False)
print(out.pivot_table(index="group",columns="gene",values="pct").to_string())
print(f"\nwrote {OUT}")
