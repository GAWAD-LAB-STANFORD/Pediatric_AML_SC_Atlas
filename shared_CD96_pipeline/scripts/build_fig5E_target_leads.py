#!/usr/bin/env python
# ====================================================================
# build_fig5E_target_leads.py  |  Figure 5 (builder)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_fig5E_target_leads.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""ALTERNATE Figure 5 panel E — TARGET bulk expression for panel-E genes:
   4 leads (CD96, CD9, SUCNR1, IL1RAP) + FLT3 (control) + 3 clinical (CD33, CD123, CLL-1).
Pulls per-sample counts from the PARENT TARGET STAR matrix (Ensembl-indexed), restricted to
the same >=80%-blast cohort the published panel used, computes log2(CPM+1), validates CD96.
Overwrites: source_data/fig5F_violin.csv, source_data/fig5F_meta_all.csv (with a `group` column).
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
PARENT = os.path.abspath(os.path.join(ROOT, "..", "ALSF_AML_plot", "H5AD"))
SD = os.path.join(ROOT, "source_data")
STAR_F = os.path.join(PARENT, "NCBI_TARGET_Star_counts_filtered.csv")

# symbol -> (ensembl, display label, group). ENSG verified vs CZ Census resolver.
GENES = [
    ("CD96",    "ENSG00000153283", "CD96",   "lead"),
    ("CD9",     "ENSG00000010278", "CD9",    "lead"),
    ("SUCNR1",  "ENSG00000198829", "SUCNR1", "lead"),
    ("IL1RAP",  "ENSG00000196083", "IL1RAP", "lead"),
    ("FLT3",    "ENSG00000122025", "FLT3",   "control"),
    ("CD33",    "ENSG00000105383", "CD33",   "clinical"),
    ("IL3RA",   "ENSG00000185291", "CD123",  "clinical"),
    ("CLEC12A", "ENSG00000172322", "CLL-1",  "clinical"),
    # PPAC-specific markers (S17 group; purple) — mark individual PPAC subsets, not pan-AML
    ("CD7",     "ENSG00000173762", "CD7",     "ppac"),
    ("TNFRSF4", "ENSG00000186827", "TNFRSF4", "ppac"),
    ("ABCA7",   "ENSG00000064687", "ABCA7",   "ppac"),
    ("ITGAX",   "ENSG00000140678", "ITGAX",   "ppac"),
]
# scRNA normal-marrow toxicity (max %HSPC, %Myeloid positive) from the count matrix
HSPC_MYE = {"CD96": 3.02, "CD9": 7.61, "SUCNR1": 7.12, "IL1RAP": 9.65,
            "FLT3": 53.06, "CD33": 28.71, "IL3RA": 12.62, "CLEC12A": 29.67,
            "CD7": 9.09, "TNFRSF4": 0.56, "ABCA7": 9.77, "ITGAX": 6.63}
TARGETABLE = 5.0

cohort = list(pd.unique(pd.read_csv(os.path.join(SD, "fig5F_violin.csv"))["sample"]))
print(f"cohort: {len(cohort)} samples")
star = pd.read_csv(STAR_F, index_col=0, low_memory=False)
have = [s for s in cohort if s in star.columns]
lib = star.loc[:, have].apply(pd.to_numeric, errors="coerce").sum(0).to_numpy(dtype=float)

absent = [g for g, e, *_ in GENES if e not in star.index]
if absent:
    print(f"  WARNING absent from STAR: {absent}")

rows, meta = [], []
for sym, ens, lab, grp in GENES:
    if ens not in star.index:
        meta.append(dict(gene=sym, label=lab, group=grp, hspc_mye=HSPC_MYE[sym],
                         toxic=bool(HSPC_MYE[sym] > 5), pct_targetable=np.nan, median=np.nan, n=0))
        continue
    counts = pd.to_numeric(star.loc[ens, have], errors="coerce").to_numpy(dtype=float)
    v = np.log2(counts / lib * 1e6 + 1.0)
    for s, val in zip(have, v):
        rows.append(dict(sample=s, gene=sym, log2cpm=float(val), label=lab, group=grp))
    meta.append(dict(gene=sym, label=lab, group=grp, hspc_mye=HSPC_MYE[sym],
                     toxic=bool(HSPC_MYE[sym] > 5),
                     pct_targetable=round(100 * float((v > TARGETABLE).mean()), 2),
                     median=round(float(np.median(v)), 4), n=int(len(v))))

pd.DataFrame(rows).to_csv(os.path.join(SD, "fig5F_violin.csv"), index=False)
m = pd.DataFrame(meta)
m.to_csv(os.path.join(SD, "fig5F_meta_all.csv"), index=False)
print(m.to_string(index=False))
cd = m.loc[m.gene == "CD96"].iloc[0]
print(f"\n  CD96 validation vs published (median 4.94, %targetable 49.9): "
      f"median {cd['median']:.2f}, %targetable {cd['pct_targetable']:.1f}")
