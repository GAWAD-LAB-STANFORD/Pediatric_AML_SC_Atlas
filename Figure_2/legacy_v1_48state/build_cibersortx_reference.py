#!/usr/bin/env python
"""Build the CIBERSORTx single-cell reference for the 30-state deconvolution.
Classes = LS_1..LS_30 (res-4.5 stable, core cells) + 0_HSPC + Myeloid_Pro (normal
refs, matching the proven Job4 run). <=300 cells/class, HVG genes, CPM. Output =
tab-delimited .txt with 'Gene' top-left and class labels repeated across the header
(the format CIBERSORTx 'Create Signature Matrix' expects)."""
import numpy as np, pandas as pd, h5py
from scipy.sparse import csr_matrix
SC="/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
OUT="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/cibersortx_input"
H5="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
CAP=300; rng=np.random.default_rng(0)
f=h5py.File(H5,"r")
def col(n):
    nd=f["obs"][n]
    if isinstance(nd,h5py.Group):
        cats=[c.decode() if isinstance(c,bytes) else c for c in nd["categories"][:]]; return np.array([cats[i] if i>=0 else "NA" for i in nd["codes"][:]],dtype=object)
    a=nd[:]; return np.array([x.decode() for x in a],dtype=object) if a.dtype.kind=="S" else a
sub=col("AML Sub-Clusters"); ct=col("Cell Type"); ci=col("clusters_for_interaction")
genes=np.array([g.decode() if isinstance(g,bytes) else g for g in f["var"]["_index"][:]])
hv=f["var"]["highly_variable"][:]; hvidx=np.where(hv)[0]; hvgenes=genes[hvidx]
N=len(sub); mask=sub!="NA"; maskidx=np.where(mask)[0]
ref=pd.read_csv(f"{SC}/jaccard_hi/reference_assignments.csv")["res4.50"].astype(int).values
pc=pd.read_csv(f"{SC}/jaccard_hi/jaccard_percluster.csv"); pc=pc[pc.res==4.5]
stab=pc[pc.mean_jaccard>=0.35].sort_values("n_cells",ascending=False); raws=stab.cluster.astype(int).tolist()
# global label array: LS states (leukemic core) + mature normal references (non-leukemic)
lab=np.array([None]*N,dtype=object)
for i,rw in enumerate(raws):
    gid=maskidx[ref==rw]; lab[gid]=f"LS_{i+1}"
nonleuk=~mask
NORM={
 "NORM_Monocyte":["CD14_Monocyte","CD16_Monocyte","AML-CD14"],
 "NORM_T_NK":["Naïve_CD4T","Naïve_CD8T","CTL","NK","Activated_CD4T","AML-CD4T","AML-CTL","AML-NK","AML-Naïve_CD8T"],
 "NORM_B":["CD20+B","PlasmaB","ProB","PreB","AML-B"],
 "NORM_Erythroid":["Erythrocytes","AML-Ery"],
 "NORM_DC":["mDC","pDC","AML-CD1C"],
 "NORM_HSPC":["Myeloid_Pro"]}
for name,cts in NORM.items():
    m=nonleuk & np.isin(ct,cts) & (lab==None)
    lab[m]=name
classes=[f"LS_{i+1}" for i in range(len(raws))]+list(NORM.keys())
# sample <=CAP per class, grouped order
sel=[]; labels=[]
for c in classes:
    idx=np.where(lab==c)[0]
    if len(idx)>CAP: idx=rng.choice(idx,CAP,replace=False)
    sel.append(idx); labels+= [c]*len(idx)
sel=np.concatenate(sel); pos={int(g):k for k,g in enumerate(sel)}
print(f"classes {len(classes)}, cells {len(sel)}, genes {len(hvidx)}",flush=True)
# stream CPM for HVG on selected cells
RC=f["layers"]["Raw_Counts"]; dd,ii,ip=RC["data"],RC["indices"],RC["indptr"][:]; G=len(genes)
mat=np.zeros((len(sel),len(hvidx)),dtype=np.float32)
selset=set(sel.tolist())
for a in range(0,N,8000):
    b=min(a+8000,N)
    chunk_sel=[g for g in range(a,b) if g in selset]
    if not chunk_sel: continue
    s,e=int(ip[a]),int(ip[b])
    X=csr_matrix((dd[s:e].astype(np.float32),ii[s:e],(ip[a:b+1]-ip[a]).astype(np.int64)),shape=(b-a,G))
    rs=np.asarray(X.sum(1)).ravel(); rs[rs==0]=1
    loc=[g-a for g in chunk_sel]
    cpm=X[loc][:,hvidx].multiply((1e6/rs[loc])[:,None]).toarray()
    for k,g in enumerate(chunk_sel): mat[pos[g]]=cpm[k]
f.close()
# write reference: genes x cells, header 'Gene' + labels
df=pd.DataFrame(mat.T, index=hvgenes, columns=labels)
path=f"{OUT}/ALSF_AML_48state_reference_CIBERSORTx.txt"
df.to_csv(path, sep="\t", index_label="Gene", float_format="%.4g")
import os; print(f"wrote {path}  ({os.path.getsize(path)/1e6:.0f} MB)")
print("class cell counts:", pd.Series(labels).value_counts().to_dict())
