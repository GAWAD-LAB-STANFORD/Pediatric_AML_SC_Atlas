#!/usr/bin/env python
"""48-state scheme (res 4.5, Jaccard >= 0.35). Label by differentiation program
ONLY (no cytogenetics), and quantify cross-cytogenetic breadth per state (the
point: shared transcriptional states that span multiple genetic subtypes)."""
import numpy as np, pandas as pd, h5py
SC="/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
OUT=f"{SC}/ls48"; import os; os.makedirs(OUT,exist_ok=True)
genes=np.load(f"{SC}/signatures_all/genes.npy",allow_pickle=True).astype(str)
M=np.load(f"{SC}/signatures_all/meancpm_res4.50.npy")            # 63 x G
jac=pd.read_csv(f"{SC}/jaccard_hi/jaccard_percluster.csv"); jac=jac[jac.res==4.5]
stab=jac[jac.mean_jaccard>=0.35].sort_values("n_cells",ascending=False)
raws=stab.cluster.astype(int).tolist(); K=len(raws)
ls=[f"LS_{i+1}" for i in range(K)]; raw2ls={rw:i for i,rw in enumerate(raws)}
print(f"{K} states at Jaccard>=0.35")

# cytogenetics + patient per cell (mask order)
f=h5py.File("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad","r")
def col(n):
    nd=f["obs"][n]; cats=[c.decode() if isinstance(c,bytes) else c for c in nd["categories"][:]]
    return np.array([cats[i] if i>=0 else "NA" for i in nd["codes"][:]],dtype=object)
sub=col("AML Sub-Clusters"); mask=sub!="NA"; cyto=col("Cytogenetic")[mask]; sid=col("SampleID")[mask]; f.close()
raw=pd.read_csv(f"{SC}/jaccard_hi/reference_assignments.csv")["res4.50"].astype(int).values

# lineage programs (curated), z-scored across the 48
PROG={"HSC/MPP":["CD34","HLF","AVP","CRHBP","HOXA9","MEIS1","MLLT3","PROM1","MECOM","GATA2"],
 "GMP":["MPO","ELANE","AZU1","PRTN3","CTSG","CEBPE","LYST","CFD"],
 "Monocyte":["LYZ","CD14","CSF1R","FCN1","VCAN","S100A8","S100A9","CD68"],
 "Macrophage":["C1QA","C1QB","C1QC","MRC1","SIGLEC1","TREM1"],
 "cDC":["CD1C","FCER1A","CLEC10A","CD207","CD1E","IRF8"],"pDC":["IRF7","LILRA4","GZMB","IL3RA"],
 "Megakaryocyte":["PPBP","PF4","ITGA2B","GP9","GP1BA","VWF","TUBB1","PF4V1"],
 "Erythroid":["GATA1","KLF1","HBB","HBA1","GYPA","ALAS2","TFRC","AHSP"],
 "Neutrophil":["FUT4","CEACAM8","MMP8","CD177","LTF","LCN2","CAMP","MMP9"],
 "B/plasma":["CD19","MS4A1","CD79A","CD79B","IGHM","JCHAIN","BLNK"],
 "T/NK":["CD3D","CD3E","IL7R","NKG7","GNLY","CD2","GZMH"],
 "Cell-cycle":["MKI67","TOP2A","CDK1","UBE2C","CENPA","BIRC5","PCNA","NUSAP1","BUB1B"]}
gi={g:i for i,g in enumerate(genes)}
Ms=np.log1p(M[raws]); Z=(Ms-Ms.mean(0))/(Ms.std(0)+1e-9)
S=pd.DataFrame({p:Z[:,[gi[g] for g in gs if g in gi]].mean(1) for p,gs in PROG.items()},index=ls)
top=S.idxmax(1); tv=S.max(1); sv=S.apply(lambda r:r.sort_values().iloc[-2],axis=1)
conf=(tv>0.6)&((tv-sv)>0.25)
lineage=np.where(conf,top,"unresolved")

rows=[]
for i,rw in enumerate(raws):
    x=ls[i]; m=raw==rw; n=m.sum()
    cyc=pd.Series(cyto[m]); cyc=cyc[~cyc.str.contains("HealthyBM")]      # leukemic cytogenetics only
    vc=cyc.value_counts(normalize=True)
    n_sub=int((vc>=0.10).sum())                                          # subtypes >=10% of state
    eff=float(np.exp(-(vc*np.log(vc)).sum())) if len(vc) else 0.0        # effective # subtypes
    topcy=f"{vc.index[0]}({vc.iloc[0]:.0%})" if len(vc) else "-"
    ps=pd.Series(sid[m]).value_counts(normalize=True); npat=(ps>=0.05).sum(); toppf=ps.iloc[0]
    pclass="single-pt" if toppf>=0.8 else ("multi-pt" if npat>=3 else "few-pt")
    apl=(cyto[m]=="PML/RARA").mean()
    rows.append([x,raw,int(n),lineage[i],round(float(tv.iloc[i]),2),n_sub,round(eff,1),topcy,pclass,int(npat),round(apl,2)])
t=pd.DataFrame(rows,columns=["LS","_raw","n_cells","lineage_program","prog_score","n_cyto>=10%","eff_n_cyto","top_cyto","patient_class","n_patients","apl_frac"])
t.drop(columns="_raw").to_csv(f"{OUT}/LS48_annotation.csv",index=False)
pd.set_option("display.width",220)
print(t.drop(columns=["_raw","prog_score","top_cyto"]).to_string(index=False))
print(f"\nlineage programs (confident):", t[t.lineage_program!='unresolved'].lineage_program.value_counts().to_dict())
print(f"cross-cytogenetic states (>=3 subtypes at >=10%):", int((t['n_cyto>=10%']>=3).sum()), "of", K)
aplst=t.sort_values('apl_frac',ascending=False).iloc[0]
print(f"APL recovered: {aplst.LS} (lineage {aplst.lineage_program}) is {aplst.apl_frac:.0%} APL, spanning {aplst['n_cyto>=10%']} subtypes >=10%")
PY
