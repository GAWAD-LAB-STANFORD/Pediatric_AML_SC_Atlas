#!/usr/bin/env python
# ====================================================================
# build_ppac_headtohead.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_ppac_headtohead.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Per-PPAC head-to-head (like Figure 5C), six panels.

For each group G (G = PPAC_1..5, or PPAC_1-5 pooled) plot every marker as a single
point in the SAME toxicity-vs-efficacy space as Figure 5C, but with efficacy made
PPAC-specific:
    x = leak     = % of normal HSPC / Myeloid_Pro positive (max of the two)   [toxicity]
    y = cov_pct  = % of G's cells positive                                    [efficacy, PPAC-specific]
  size = rec_pct = % of contributing patients (>=20 cells in G) with >=20% of
                   their G cells positive                                     [recurrence]

Markers per panel = the TOP-5 PPAC-specific markers for G, ranked by % of the
PPAC's cells positive (from funnel_ppac_specific.csv, the S17 ranking; toxicity
tolerance <=10% on both normal compartments) UNION all 9 clinical-stage AML targets.
Clinical targets keep their fixed identity colour; discovered markers are grey.
This shows, per PPAC, where the newly-discovered surface markers AND the clinical
targets fall on the efficacy/toxicity plane.

Output -> source_data/figS_ppac_headtohead.csv  (long: group x marker).
Builder reads the bundled count matrix only.
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
CM = os.path.join(HERE, "..", "count_matrix", "ALSF_AML_raw_counts.h5ad")
RANK = os.path.join(HERE, "..", "source_data", "funnel_ppac_specific.csv")
COMP = os.path.join(HERE, "..", "source_data", "funnel_ppac_specific_composition.csv")
OUT = os.path.join(HERE, "..", "source_data", "figS_ppac_headtohead.csv")

REC_THRESH = 0.20
MIN_CELLS_PER_PT = 20
TOP_K = 5
PPACS = ["PPAC_1", "PPAC_2", "PPAC_3", "PPAC_4", "PPAC_5"]
# the 12 panel-E antigens (= the genes in Figure 7B/C): leads + FLT3 control + clinical + PPAC-specific
CLIN = {"CD96": "CD96", "CD9": "CD9", "SUCNR1": "SUCNR1", "IL1RAP": "IL1RAP", "FLT3": "FLT3",
        "CD33": "CD33", "IL3RA": "CD123", "CLEC12A": "CLL-1", "CD7": "CD7", "TNFRSF4": "TNFRSF4",
        "ABCA7": "ABCA7", "ITGAX": "ITGAX"}

# ---- pick top-5 discovered markers per group from the S17 ranking (% positive) ----
rank = pd.read_csv(RANK)
top_by_group = {}
for g, sub in rank.groupby("group"):
    top_by_group[g] = list(sub.sort_values("cov_pct", ascending=False)
                           .head(TOP_K)["gene"])
GROUPS = ["All PPAC"] + PPACS

# union of every gene we need a column for
need = set(CLIN)
for g in GROUPS:
    need |= set(top_by_group.get(g, []))
need = sorted(need)

print("Loading count matrix ...")
A = ad.read_h5ad(CM)
present = [g for g in need if g in A.var_names]
missing = [g for g in need if g not in A.var_names]
if missing:
    print(f"  WARNING: not in matrix (skipped): {missing}")
pac = A.obs["Prognosis-Associated Clusters"].astype(str).values
samp = A.obs["Sample"].astype(str).values
hspc = pac == "0_HSPC"
mye = pac == "Myeloid_Pro"
X = A[:, present].X
pos = (np.asarray(X.todense()) > 0) if sp.issparse(X) else (np.asarray(X) > 0)
present = np.array(present)
col = {g: i for i, g in enumerate(present)}

aml_samples = sorted({s for s in np.unique(samp) if "Healthy" not in s})
hspc_pct = 100 * pos[hspc].mean(0)
mye_pct = 100 * pos[mye].mean(0)
comp = pd.read_csv(COMP).set_index("group")

rows = []
for gname in GROUPS:
    gmask = (pac == gname) if gname != "All PPAC" else np.isin(pac, PPACS)
    n_cells = int(gmask.sum())
    contrib = [s for s in aml_samples if (gmask & (samp == s)).sum() >= MIN_CELLS_PER_PT]
    spflag = bool(comp.loc[gname, "single_patient_flag"]) if gname in comp.index else False
    markers = list(dict.fromkeys(top_by_group.get(gname, []) + list(CLIN)))  # top5 then clinical, dedup
    cov_all = 100 * pos[gmask].mean(0)
    for g in markers:
        if g not in col:
            continue
        i = col[g]
        cov = float(cov_all[i])
        if contrib:
            rec = 100 * np.mean([pos[gmask & (samp == s)][:, i].mean() >= REC_THRESH
                                 for s in contrib])
        else:
            rec = 0.0
        leak = float(max(hspc_pct[i], mye_pct[i]))
        rows.append(dict(
            group=gname, gene=g, label=CLIN.get(g, g),
            cls="clinical" if g in CLIN else "discovered",
            leak=round(leak, 2), cov_pct=round(cov, 2), rec_pct=round(float(rec), 2),
            HSPC=round(float(hspc_pct[i]), 2), Myeloid=round(float(mye_pct[i]), 2),
            in_top5=g in top_by_group.get(gname, []),
            n_cells=n_cells, single_patient_flag=spflag))

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(f"  wrote {OUT}  ({len(df)} rows, {df.group.nunique()} panels)\n")

pd.set_option("display.width", 170)
for gname in GROUPS:
    sub = df[df.group == gname].sort_values(["cls", "cov_pct"], ascending=[True, False])
    tag = "  [single-patient PPAC]" if sub.single_patient_flag.iloc[0] else ""
    print(f"=== {gname}{tag} ===")
    print(sub[["label", "cls", "in_top5", "leak", "cov_pct", "rec_pct"]].to_string(index=False))
    print()
