#!/usr/bin/env python
"""Stitch the 4 CIBERSORTx 30-state batches (Job14-17), renormalize LS fractions
among leukemic states, join TARGET survival + cytogenetics, and run cross-patient
detection + cytogenetic sanity. Writes survival input for the Cox step."""
import numpy as np, pandas as pd
PKG="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated"; OUT=f"{PKG}/cibersortx_output"
SL="/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
d=pd.concat([pd.read_csv(f"{OUT}/CIBERSORTx_Job{j}_Results.csv") for j in (14,15,16,17)],ignore_index=True)
print(f"stitched {len(d)} patients; median correlation {d['Correlation'].median():.3f}, median RMSE {d['RMSE'].median():.3f}")
ls=[f"LS_{i}" for i in range(1,31)]; norm=[c for c in d.columns if c.startswith("NORM_")]
d["leukemic_purity"]=d[ls].sum(1)/(d[ls].sum(1)+d[norm].sum(1)).replace(0,1)
frac=d[["Mixture"]+ls].copy(); frac.columns=["sample"]+ls
lss=frac[ls].values; s=lss.sum(1,keepdims=True); s[s==0]=1
frac[ls]=lss/s   # renormalize among leukemic states
frac["leukemic_purity"]=d["leukemic_purity"].values
frac.to_csv(f"{SL}/target_ls/LS_target_fractions_cibersortx.csv",index=False)
surv=pd.read_csv(f"{PKG}/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/../../08_Source_Data/Figure_2/Figure2E_TARGET_survival_input.csv") if False else pd.read_csv(f"{PKG}/__SUBMISSION_PACKAGE/08_Source_Data/Figure_2/Figure2E_TARGET_survival_input.csv")
cd=pd.read_csv(f"{PKG}/__SUBMISSION_PACKAGE/08_Source_Data/Figure_2/Figure2CD_TARGET_PAC_perpatient.csv")[["sample","subtype2"]]
J=frac.merge(surv[["sample","os_time","os_event","efs_time","efs_event","age","wbc","risk"]],on="sample",how="left").merge(cd.drop_duplicates("sample"),on="sample",how="left")
J.to_csv(f"{SL}/target_ls/LS_target_survival_cibersortx.csv",index=False)
print(f"mean leukemic purity {frac['leukemic_purity'].mean():.2f}; patients with survival {J.os_time.notna().sum()}, deaths {int(J.os_event.sum())}")
# cross-patient detection + cytogenetic sanity
disc=pd.read_csv(f"{SL}/ls_scheme/LS_discovery_validation.csv").set_index("LS")
alias={"CBFB/MYH11":"inv(16)","RUNX1/RUNX1T1":"t(8;21)","MLLr":"KMT2A(MLL)","NUP98/NSD1":"NUP98-r","PML/RARA":"PML-RARA"}
jc=J.dropna(subset=["subtype2"]); cyto=jc.groupby("subtype2")[ls].mean().T
ok=0; det=[]
print("\n=== cytogenetic-identity sanity (discovery -> TARGET peak) ===")
for x in ls:
    dcy=str(disc.loc[x,"top_cyto(enrich)"]).split("(")[0] if x in disc.index else "-"
    dexp=alias.get(dcy,dcy); peak=cyto.loc[x].idxmax(); m=(dexp==peak)and dcy not in ("-","nan"); ok+=m
    det.append([x,round((frac[x]>0.01).mean(),3),round((frac[x]>0.05).mean(),3),round(float(frac[x].mean()),4)])
    print(f"  {x}: disc[{dcy}] peak[{peak} {cyto.loc[x].max():.3f}] {'OK' if m else ''}")
print(f"\nSANITY: {ok}/30 states recover their discovery cytogenetic identity in TARGET")
pd.DataFrame(det,columns=["LS","det>1%","det>5%","mean_frac"]).to_csv(f"{SL}/target_ls/LS_cibersortx_detection.csv",index=False)
