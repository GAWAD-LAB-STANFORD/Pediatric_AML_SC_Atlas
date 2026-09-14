#!/usr/bin/env python
"""S24 panels E/F data — CD96–nectin checkpoint axis (both directions).
CD96 is aberrant on AML (target) AND is the inhibitory receptor on T/NK. Builds, for
AML (4-type leukemic) + the 12 T/NK subsets: % positive + mean log1p of the PVR-axis
receptors (CD96/TIGIT/CD226) and nectin ligands (PVR/NECTIN1/NECTIN2), then the
directional co-expression connectivity (% ligand+ sender × % CD96+ receiver ÷ 100) for
AML→T/NK vs T/NK→AML. Atlas + T/NK object; raw counts; no fabricated values.
Outputs: source_data/cd96_axis_expr.csv, source_data/cd96_axis_connectivity.csv
"""
import h5py, numpy as np, pandas as pd, os
from scipy.sparse import csr_matrix
dec=lambda a:np.array([x.decode() if isinstance(x,bytes) else x for x in a])
ROOT="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated"
RD=os.path.join(ROOT,"__SUBMISSION_PACKAGE_v2_recluster/07_Code/shared_CD96_pipeline/source_data")
H5=os.path.join(ROOT,"ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad")
H5T=os.path.join(ROOT,"ALSF_AML_plot/H5AD/ALSF_AML_total_Tcell_rna_withTCR.h5ad")
PC=os.path.join(ROOT,"__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2/tnk_subset_percell.csv")
GENES=["CD96","TIGIT","CD226","PVR","NECTIN1","NECTIN2"]
ROLE={"CD96":"receptor","TIGIT":"receptor","CD226":"receptor","PVR":"ligand","NECTIN1":"ligand","NECTIN2":"ligand"}
def raw(H5):
    f=h5py.File(H5,"r"); var=dec(f["raw"]["var"]["_index"][:]); gp={g:i for i,g in enumerate(var)}
    gi=[gp[g] for g in GENES if g in gp]; grp=f["raw"]["X"]
    X=csr_matrix((grp["data"][:],grp["indices"][:],grp["indptr"][:]),shape=tuple(grp.attrs["shape"]) if "shape" in grp.attrs else None)
    sub=np.asarray(X[:,gi].todense()); f.close(); return sub,[g for g in GENES if g in gp]
rows=[]
f=h5py.File(H5,"r"); obs=dec(f["obs"]["_index"][:]); g=f["obs"]["Cell Type"]; ct=dec(g["categories"][:])[g["codes"][:]]; f.close()
Xa,pa=raw(H5); LEUK=np.isin(ct,["AML","AML-MKI67","AML-PCNA","AML-CD1C"])
for i,gene in enumerate(pa):
    v=Xa[LEUK,i]; rows.append(dict(compartment="AML (leukemic)",gene=gene,role=ROLE[gene],pct=round(100*(v>0).mean(),1),mean=round(np.log1p(v).mean(),3)))
ft=h5py.File(H5T,"r"); bc=dec(ft["obs"]["_index"][:]); ft.close(); Xt,pt=raw(H5T)
pc=pd.read_csv(PC); idx=pd.Series(np.arange(len(bc)),index=bc)
ORD=["Naïve CD4 T","Memory CD4 T","Treg","Naïve CD8 T","Memory CD8 T","GZMK CD8 T","GZMB CD8 T","GZMB DNT","MAIT","GZMK NK","GZMB NK","Proliferating T"]
for s in ORD:
    ii=idx.reindex(pc.loc[pc.subset==s,"barcode"]).dropna().astype(int).values
    for i,gene in enumerate(pt):
        v=Xt[ii,i]; rows.append(dict(compartment=s,gene=gene,role=ROLE[gene],pct=round(100*(v>0).mean(),1),mean=round(np.log1p(v).mean(),3)))
d=pd.DataFrame(rows); d.to_csv(os.path.join(RD,"cd96_axis_expr.csv"),index=False)
P=d.pivot_table(index="compartment",columns="gene",values="pct"); aml=P.loc["AML (leukemic)"]
crows=[]
for s in ORD:
    t=P.loc[s]
    crows.append(dict(subset=s,direction="AML → T/NK (nectin → CD96)",score=round((aml["PVR"]+aml["NECTIN1"])*t["CD96"]/100,2)))
    crows.append(dict(subset=s,direction="T/NK → AML (nectin → CD96)",score=round((t["PVR"]+t["NECTIN1"])*aml["CD96"]/100,2)))
pd.DataFrame(crows).to_csv(os.path.join(RD,"cd96_axis_connectivity.csv"),index=False)
print("wrote cd96_axis_expr.csv + cd96_axis_connectivity.csv")
