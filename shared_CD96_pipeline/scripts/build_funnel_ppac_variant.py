#!/usr/bin/env python
# ====================================================================
# build_funnel_ppac_variant.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_funnel_ppac_variant.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Discovery-funnel variant: step 4 = 'any PPAC > 20%' instead of 'AML > 20% in >=1 patient'.
Reproduces the surfaceome -> HSPC<3% -> Myeloid<=5% sparing gates from the count matrix, then
applies the two alternative step-4 criteria and reports what comes back. Output ->
source_data/funnel_ppac_variant.csv + printed summary."""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
CM = os.path.join(HERE, "..", "count_matrix", "ALSF_AML_raw_counts.h5ad")
SURF = os.path.join(HERE, "..", "..", "ref", "surfaceome_genes.txt")
OUT = os.path.join(HERE, "..", "source_data", "funnel_ppac_variant.csv")
GATE_HSPC, GATE_MYE, PPAC_MIN, AML_MIN = 3.0, 5.0, 20.0, 20.0
PPACS = ["PPAC_1", "PPAC_2", "PPAC_3", "PPAC_4", "PPAC_5"]
LEUK = [f"AML_{i}" for i in range(1, 28)]
CLIN = {"CD96", "CD33", "IL3RA", "CLEC12A", "FLT3", "MSLN", "CD70", "ADGRG1", "IL1RAP"}

A = ad.read_h5ad(CM)
surf = [g.strip() for g in open(SURF) if g.strip()]
present = [g for g in surf if g in A.var_names]
pac = A.obs["Prognosis-Associated Clusters"].astype(str).values
sub = A.obs["AML Sub-Clusters"].astype(str).values
samp = A.obs["Sample"].astype(str).values
leuk = np.isin(sub, LEUK)
hspc = pac == "0_HSPC"; mye = pac == "Myeloid_Pro"
X = A[:, present].X
pos = (np.asarray(X.todense()) > 0) if sp.issparse(X) else (np.asarray(X) > 0)   # cells x Npresent

hspc_pct = 100 * pos[hspc].mean(0)
mye_pct = 100 * pos[mye].mean(0)
aml_pct = 100 * pos[leuk].mean(0)
ppac_pct = {p: 100 * pos[pac == p].mean(0) for p in PPACS}
maxppac = np.max([ppac_pct[p] for p in PPACS], axis=0)
which = np.array(PPACS)[np.argmax([ppac_pct[p] for p in PPACS], axis=0)]
# AML > 20% in >=1 patient (the current step-4 metric)
samples = [s for s in np.unique(samp) if "Healthy" not in s]
npts = np.zeros(len(present), int)
for s in samples:
    m = leuk & (samp == s)
    if m.sum():
        npts += (pos[m].mean(0) >= AML_MIN / 100)
df = pd.DataFrame(dict(gene=present, HSPC=hspc_pct.round(2), Myeloid=mye_pct.round(2),
                       AML=aml_pct.round(2), max_PPAC=maxppac.round(2), max_PPAC_cluster=which,
                       n_pts_AML20=npts, clinical=[g in CLIN for g in present]))
df.to_csv(OUT, index=False)

spar = df[(df.HSPC < GATE_HSPC) & (df.Myeloid <= GATE_MYE)]
old4 = spar[spar.n_pts_AML20 >= 1]                       # current step 4
new4 = spar[spar.max_PPAC > PPAC_MIN]                     # requested step 4

print("=== FUNNEL ===")
print(f"  surfaceome detectable in scRNA : {len(df)}")
print(f"  HSPC-sparing (<3%)             : {(df.HSPC<GATE_HSPC).sum()}")
print(f"  + Myeloid-sparing (<=5%)       : {len(spar)}")
print(f"  step4 [AML>20% in >=1 pt]      : {len(old4)}   (current funnel)")
print(f"  step4 [ANY PPAC>20%]           : {len(new4)}   (requested)")
print(f"  overlap of the two step-4 sets : {len(set(old4.gene)&set(new4.gene))}")
print(f"\n=== genes from the NEW step4 (any PPAC>20%), top 40 by max_PPAC ===")
show = new4.sort_values("max_PPAC", ascending=False)
print(show[["gene","HSPC","Myeloid","AML","max_PPAC","max_PPAC_cluster","n_pts_AML20","clinical"]].head(40).to_string(index=False))
print(f"\nclinical targets in the NEW set: {sorted(new4[new4.clinical].gene)}")
print(f"in NEW but NOT in current AML set (novel via PPAC): {sorted(set(new4.gene)-set(old4.gene))[:40]}")
print(f"in current AML set but dropped by PPAC criterion : {sorted(set(old4.gene)-set(new4.gene))[:40]}")
