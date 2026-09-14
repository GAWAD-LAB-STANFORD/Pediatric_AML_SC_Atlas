#!/usr/bin/env python
"""Figure 3 step 6: build the mean-expression matrices for panel C (focus-state
identity + targetable surface antigens vs normal HSPC / favorable / poor) and
panel D (HOX expression across CBF-prognostic states). Uses per-state mean CPM +
a normal-HSPC reference (0_HSPC) from the h5ad."""
import os, sys, re
import scanpy as sc, h5py, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
genes=np.load(C.GENES,allow_pickle=True).astype(str); M=np.load(C.MEANCPM)
raws,ls=C.stable_raws(); Ms=pd.DataFrame(np.log1p(M[raws]),index=ls,columns=genes)
order=pd.read_csv(C.CBF_ORDER); focus=open(os.path.join(C.DATA,"LS_focus_state.txt")).read().strip()
cbf_states=order.LS.tolist()
cox=pd.read_csv(os.path.join(C.DATA,"LS_cbf_cox_efs.csv"))
favs=order.LS[order.dir=="favorable"].tolist(); poors=[x for x in order.LS[order.dir=="poor"].tolist() if x!=focus]
# normal HSPC reference
f=h5py.File(C.H5,"r"); nd=f["obs"]["clusters_for_interaction"]
cats=[c.decode() if isinstance(c,bytes) else c for c in nd["categories"][:]]
ci=np.array([cats[i] if i>=0 else "NA" for i in nd["codes"][:]],dtype=object); f.close()
# focus-state signature = its OWN top up-genes (data-driven), not a hardcoded set
up=pd.read_csv(os.path.join(C.DATA,"LS_upgenes_wilcoxon.csv"))
fsig=[g for g in up[up.LS==focus].sort_values("log2FC",ascending=False).gene.tolist() if g in Ms.columns][:8]
groups={"not LT-HSC":["HLF","AVP","SPINK2"],
 "focus signature":fsig,
 "shared progenitor":["CD34","GATA2","MEIS1","HOXA9"],
 "not GMP (granule)":["MPO","PRTN3","ELANE","AZU1"],"cycling":["MKI67","TOP2A","CDK1"],
 "surface targets":["CD96","FLT3","CD38","IL3RA","CLEC12A","CD33","KIT","CRLF2","ITGAX","TNFRSF4"]}
mk=[g for gs in groups.values() for g in gs]
ad=sc.read_h5ad(C.H5,backed="r")[np.where(ci=="0_HSPC")[0]].to_memory()
lyr="Raw_Counts" if "Raw_Counts" in ad.layers else "counts"; ad.X=ad.layers[lyr].copy(); sc.pp.normalize_total(ad,target_sum=1e6)
present=[g for g in mk if g in ad.var_names and g in Ms.columns]
X=ad[:,present].X; X=X.toarray() if hasattr(X,"toarray") else np.asarray(X); nh=pd.Series(np.log1p(X.mean(0)),index=present)
pc=pd.DataFrame({"gene":present,"group":[next(grp for grp,gs in groups.items() if g in gs) for g in present],
  "Normal_HSPC":nh.values,"Favorable_mean":[Ms.loc[favs,g].mean() if favs else np.nan for g in present],
  "Focus":[Ms.loc[focus,g] for g in present],"Poor_mean":[Ms.loc[poors,g].mean() if poors else np.nan for g in present]})
pc.to_csv(os.path.join(C.DATA,"panelC_targetable.csv"),index=False)
# panel D: HOX expression across CBF-prognostic states
hox=[g for g in genes if re.match(r"^HOX[A-D]\d+$",g)]
hox=sorted(hox,key=lambda g:(re.match(r"^HOX([A-D])(\d+)$",g).group(1),int(re.match(r"^HOX([A-D])(\d+)$",g).group(2))))
Ms.loc[cbf_states,hox].T.to_csv(os.path.join(C.DATA,"panelD_hox.csv"))
print(f"wrote panelC_targetable.csv ({len(present)} markers) and panelD_hox.csv ({len(hox)} HOX x {len(cbf_states)} states); focus={focus}")
