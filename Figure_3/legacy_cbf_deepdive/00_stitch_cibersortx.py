#!/usr/bin/env python
"""Figure 3 step 0: stitch the 4 CIBERSORTx job outputs (disjoint TARGET batches)
into LS_survival.csv = per-patient LS_* + NORM_* fractions (NEW deconvolution) joined
to the EXISTING, unchanged TARGET clinical annotations (subtype2, OS, EFS). No values
are altered; fractions come straight from CIBERSORTx, clinical from the TARGET mapping."""
import os, sys, glob
import pandas as pd, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
CBX = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/cibersortx_output"
JOBS = [f"{CBX}/CIBERSORTx_Job{j}_Results.csv" for j in (24,25,26,27)]
SL = "/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
CLIN = f"{SL}/target_ls/LS48_survival.csv"          # sample, os_time, os_event, subtype2 (real TARGET clinical)
EFS  = f"{C.BASE}/08_Source_Data/Figure_2/Figure2E_TARGET_survival_input.csv"  # sample, efs_time, efs_event

frames=[pd.read_csv(j) for j in JOBS]
for j,fr in zip(JOBS,frames): print(f"  {os.path.basename(j)}: {fr.shape[0]} patients, {fr.shape[1]} cols")
cbx=pd.concat(frames, ignore_index=True)
dup=cbx["Mixture"].duplicated().sum()
print(f"stitched: {cbx.shape[0]} rows, {cbx['Mixture'].nunique()} unique patients, {dup} duplicates")
assert dup==0, "duplicate patients across jobs!"
cbx=cbx.rename(columns={"Mixture":"sample"})
lscols=[c for c in cbx.columns if c.startswith("LS_")]
normcols=[c for c in cbx.columns if c.startswith("NORM_")]
print(f"leukemic-state cols: {len(lscols)}, normal cols: {len(normcols)}")
print(f"CIBERSORTx fit: median Correlation={cbx['Correlation'].median():.3f}, "
      f"patients with P-value<0.05: {(cbx['P-value']<0.05).sum()}/{len(cbx)}")
# join REAL clinical (unchanged)
clin=pd.read_csv(CLIN)[["sample","subtype2","os_time","os_event"]]
efs=pd.read_csv(EFS)[["sample","efs_time","efs_event"]]
m=cbx.merge(clin,on="sample",how="left").merge(efs,on="sample",how="left")
print(f"matched to clinical: subtype2 {m['subtype2'].notna().sum()}/{len(m)}, "
      f"OS {m['os_time'].notna().sum()}, EFS {m['efs_time'].notna().sum()}")
print(f"CBF patients (t(8;21)+inv(16)): {m['subtype2'].isin(['t(8;21)','inv(16)']).sum()}")
keep=["sample"]+lscols+normcols+["subtype2","os_time","os_event","efs_time","efs_event","P-value","Correlation","RMSE"]
out=os.path.join(C.DATA,"LS_survival.csv")
m[keep].to_csv(out,index=False)
print(f"wrote {out} ({m.shape[0]} patients x {len(lscols)} leukemic states)")
