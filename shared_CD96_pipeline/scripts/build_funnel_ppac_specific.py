#!/usr/bin/env python
# ====================================================================
# build_funnel_ppac_specific.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_funnel_ppac_specific.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""PPAC-specific discovery (supplementary, S17).

Ranks surface markers WITHIN each poor-prognosis cluster (PPAC) by the simplest
PPAC-specific efficacy metric -- the **percent of that PPAC's cells that are
positive** (no composite). Candidates first pass a toxicity funnel: detectable in
scRNA, normal HSPC <= GATE% AND normal Myeloid_Pro <= GATE%, and the TCR/CD3
family dropped (MPAL-like co-expression). Toxicity tolerance is set to 10% for
both HSPC and Myeloid_Pro (more permissive than the main Figure-5A funnel, which
uses HSPC<3% / Myeloid<=5%); this is intentional and noted in the legend.

  metric  cov_pct(gene, G) = % of group G's cells positive (count > 0)
  groups  G = PPAC_1..5, or "All PPAC" = PPAC_1-5 pooled

`rec_pct` (% of contributing patients with >=20% of their G cells positive) and a
per-PPAC patient-composition audit are also reported so single-patient PPACs are
flagged (their positivity is driven by one patient and is not generalizable).

Outputs:
  source_data/funnel_ppac_specific.csv               (long: group x gated-gene metrics)
  source_data/funnel_ppac_specific_composition.csv   (per-group patient composition)
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
CM = os.path.join(HERE, "..", "count_matrix", "ALSF_AML_raw_counts.h5ad")
SURF = os.path.join(HERE, "..", "..", "ref", "surfaceome_genes.txt")
OUT = os.path.join(HERE, "..", "source_data", "funnel_ppac_specific.csv")
OUT_COMP = os.path.join(HERE, "..", "source_data", "funnel_ppac_specific_composition.csv")

GATE_HSPC, GATE_MYE = 10.0, 10.0   # toxicity tolerance: <=10% for BOTH normal compartments
REC_THRESH = 0.20            # marker "covers" a patient's PPAC if >=20% of those cells positive
MIN_CELLS_PER_PT = 20       # patient must have >=20 cells in the group to be scored for recurrence
TOP_N = 12
PPACS = ["PPAC_1", "PPAC_2", "PPAC_3", "PPAC_4", "PPAC_5"]
# TCR / CD3 family -- MPAL-like contamination, dropped by the toxicity funnel (matches Fig 5A)
TCR_CD3_EXCLUDE = {"CD3D", "CD3E", "CD3G", "CD247",
                   "TRAC", "TRBC1", "TRBC2", "TRGC1", "TRGC2", "TRDC"}
CLIN = {"CD96", "CD33", "IL3RA", "CLEC12A", "FLT3", "MSLN", "CD70", "ADGRG1", "IL1RAP"}
GLABEL = {"IL3RA": "CD123", "CLEC12A": "CLL-1", "ADGRG1": "GPR56"}

print("Loading count matrix ...")
A = ad.read_h5ad(CM)
surf = [g.strip() for g in open(SURF) if g.strip()]
present = [g for g in surf if g in A.var_names]
pac = A.obs["Prognosis-Associated Clusters"].astype(str).values
samp = A.obs["Sample"].astype(str).values
hspc = pac == "0_HSPC"
mye = pac == "Myeloid_Pro"
X = A[:, present].X
pos = (np.asarray(X.todense()) > 0) if sp.issparse(X) else (np.asarray(X) > 0)  # cells x Npresent
present = np.array(present)
print(f"  {len(present)} surfaceome genes detectable in scRNA")

# ---- toxicity funnel (tolerance <=10% on both normal compartments) ----
hspc_pct = 100 * pos[hspc].mean(0)
mye_pct = 100 * pos[mye].mean(0)
is_tcr = np.array([g in TCR_CD3_EXCLUDE for g in present])
gate = (hspc_pct <= GATE_HSPC) & (mye_pct <= GATE_MYE) & (~is_tcr)
print(f"  HSPC-sparing (<={GATE_HSPC}%)        : {(hspc_pct <= GATE_HSPC).sum()}")
print(f"  + Myeloid-sparing (<={GATE_MYE}%)    : {((hspc_pct <= GATE_HSPC) & (mye_pct <= GATE_MYE)).sum()}")
print(f"  + drop TCR/CD3 (toxicity funnel)  : {gate.sum()}  <- genes eligible to be ranked\n")

aml_samples = sorted({s for s in np.unique(samp) if "Healthy" not in s and "HealthyBM" not in s})

groups = [(p, pac == p) for p in PPACS] + [("All PPAC", np.isin(pac, PPACS))]

comp_rows = []   # per-group composition audit
long_rows = []   # per-group x gene metrics (gated genes only)

for gname, gmask in groups:
    n_cells = int(gmask.sum())
    # patient composition of this group
    pt_counts = {s: int((gmask & (samp == s)).sum()) for s in aml_samples}
    pt_counts = {s: c for s, c in pt_counts.items() if c > 0}
    contrib = [s for s, c in pt_counts.items() if c >= MIN_CELLS_PER_PT]
    n_any = len(pt_counts)
    n_contrib = len(contrib)
    top_share = (max(pt_counts.values()) / n_cells * 100) if n_cells else 0.0
    comp_rows.append(dict(group=gname, n_cells=n_cells, n_patients_any=n_any,
                          n_patients_contrib=n_contrib,
                          dominant_patient_share_pct=round(top_share, 1),
                          single_patient_flag=(n_contrib < 3)))

    # primary metric: % of the PPAC's cells positive
    cov = 100 * pos[gmask].mean(0)
    # patient-recurrence (informational; flags single-patient signal)
    if n_contrib:
        hit = np.zeros(len(present), int)
        for s in contrib:
            m = gmask & (samp == s)
            hit += (pos[m].mean(0) >= REC_THRESH)
        rec = 100 * hit / n_contrib
    else:
        rec = np.zeros(len(present))

    for i in np.where(gate)[0]:
        g = present[i]
        long_rows.append(dict(
            group=gname, gene=g, label=GLABEL.get(g, g),
            cov_pct=round(float(cov[i]), 2),
            rec_pct=round(float(rec[i]), 2),
            HSPC=round(float(hspc_pct[i]), 2), Myeloid=round(float(mye_pct[i]), 2),
            n_contrib_patients=n_contrib,
            clinical=g in CLIN))

comp = pd.DataFrame(comp_rows)
long = pd.DataFrame(long_rows)
long.to_csv(OUT, index=False)
comp.to_csv(OUT_COMP, index=False)

pd.set_option("display.width", 160)
print("=== PER-PPAC PATIENT COMPOSITION (artifact audit) ===")
print(comp.to_string(index=False))
print("  single_patient_flag = fewer than 3 patients contribute >=20 cells "
      "-> positivity driven by one patient, not generalizable\n")

for gname, _ in groups:
    sub = long[long.group == gname].sort_values("cov_pct", ascending=False).head(TOP_N)
    flag = comp.loc[comp.group == gname, "single_patient_flag"].iloc[0]
    ncon = comp.loc[comp.group == gname, "n_patients_contrib"].iloc[0]
    tag = "  [SINGLE/FEW-PATIENT -- not generalizable]" if flag else ""
    print(f"=== {gname}  (n_contrib_patients={ncon}){tag} ===")
    show = sub[["label", "cov_pct", "rec_pct", "HSPC", "Myeloid", "clinical"]]
    print(show.to_string(index=False))
    print()
