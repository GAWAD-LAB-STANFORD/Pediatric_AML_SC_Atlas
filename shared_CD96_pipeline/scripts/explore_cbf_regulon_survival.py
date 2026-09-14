#!/usr/bin/env python
# ====================================================================
# explore_cbf_regulon_survival.py  |  Exploratory analysis (CBF regulon survival screen)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python explore_cbf_regulon_survival.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""EXPLORATION (Figure 3A idea): do the single-cell ERG / IKZF2 / LHX6 regulons stratify outcome
in CBF-AML (RUNX1-RUNX1T1 + CBFB-MYH11) in the TARGET bulk cohort?

Regulon gene sets come from the SINGLE-CELL SCENIC run (Total_30_scenic_integrated-output loom,
row_attrs['Regulons']) — the bulk SCENIC run did NOT recover these three regulons, so we score
the single-cell-derived target-gene sets on the bulk expression as a module z-score (the bulk
analog of AUCell regulon activity). Survival from Figure7F_TARGET_CD96_merged.csv.
NO fabricated values — every number derives from the loom counts + the TARGET clinical table.
-> source_data/fig3_cbf_regulon_survival.csv  + printed stats
"""
import os, h5py
import numpy as np, pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

ROOT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated"
SC_LOOM = f"{ROOT}/ALSF_AML_plot/H5AD/Total_30_scenic_integrated-output_python_3.8_new.loom"
T_LOOM  = f"{ROOT}/ALSF_AML_plot/H5AD/NCI_Target_scenic_integrated-output.loom"
CLIN    = f"{ROOT}/Figure7/Figure7F_TARGET_CD96_merged.csv"
SD      = f"{ROOT}/AML_2026_Alternate/source_data"
REGS = ["ERG", "IKZF2", "LHX6"]

# ---- 1) single-cell regulon target-gene sets ----
with h5py.File(SC_LOOM, "r") as h:
    reg = h["row_attrs"]["Regulons"]
    sc_genes = np.array([g.decode() if isinstance(g, bytes) else g for g in h["row_attrs"]["Gene"][:]])
    regsets = {}
    for tf in REGS:
        col = [n for n in reg.dtype.names if n.lower().startswith(tf.lower() + "_")][0]
        regsets[tf] = list(sc_genes[reg[col][:].astype(bool)])
for tf in REGS:
    print(f"  {tf} regulon: {len(regsets[tf])} genes -> {regsets[tf]}")

# ---- 2) TARGET bulk: CBF primary samples + expression of regulon genes ----
with h5py.File(T_LOOM, "r") as t:
    dec = lambda k: np.array([x.decode() if isinstance(x, bytes) else x for x in t["col_attrs"][k][:]])
    cid, cyto, styp = dec("CellID"), dec("Primary Cytogenetic Code"), dec("Sample.Type")
    genes = np.array([g.decode() if isinstance(g, bytes) else g for g in t["row_attrs"]["Gene"][:]])
    is_cbf = np.isin(cyto, ["RUNX1-RUNX1T1", "CBFB-MYH11"]) & np.char.startswith(styp, "Primary")
    idx = np.where(is_cbf)[0]
    sub = t["matrix"][:, idx].astype(float)                      # genes x CBF samples (raw counts)
    smp, scyto = cid[idx], cyto[idx]
print(f"\nCBF primary samples in loom: {len(idx)}  "
      f"(RUNX1-RUNX1T1 {int((scyto=='RUNX1-RUNX1T1').sum())}, CBFB-MYH11 {int((scyto=='CBFB-MYH11').sum())})")

lib = sub.sum(0); logcpm = np.log2(sub / lib * 1e6 + 1.0)        # genes x samples
g2i = {g: i for i, g in enumerate(genes)}
def module_z(geneset):
    rows = [g2i[g] for g in geneset if g in g2i]
    M = logcpm[rows, :]                                          # k genes x samples
    Z = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9)
    return Z.mean(0), len(rows), len(geneset)
df = pd.DataFrame({"sample": smp, "cyto": scyto})
for tf in REGS:
    sc, k, n = module_z(regsets[tf]); df[tf] = sc
    print(f"  {tf}: scored {k}/{n} genes present in TARGET")
df["composite"] = df[REGS].mean(1)

# ---- 3) join survival ----
cl = pd.read_csv(CLIN, low_memory=False)
keep = ["sample", "Overall Survival Time in Days", "Event Free Survival Time in Days",
        "Vital Status", "First Event"]
cl = cl[keep].drop_duplicates("sample")
df = df.merge(cl, on="sample", how="inner")
df["OS_time"] = pd.to_numeric(df["Overall Survival Time in Days"], errors="coerce") / 365.25
df["OS_event"] = (df["Vital Status"].astype(str).str.lower() == "dead").astype(int)
df["EFS_time"] = pd.to_numeric(df["Event Free Survival Time in Days"], errors="coerce") / 365.25
df["EFS_event"] = (~df["First Event"].astype(str).str.lower().isin(["censored", "none", "nan"])).astype(int)
df = df.dropna(subset=["OS_time"]).reset_index(drop=True)
df.to_csv(os.path.join(SD, "fig3_cbf_regulon_survival.csv"), index=False)
print(f"\nAnalyzable CBF patients with survival: {len(df)} "
      f"(deaths {df.OS_event.sum()}, EFS events {df.EFS_event.sum()})")

# ---- 4) survival association: Cox (continuous) + median-split log-rank ----
def analyse(d, label):
    print(f"\n===== {label}  (n={len(d)}, deaths={d.OS_event.sum()}, EFS events={d.EFS_event.sum()}) =====")
    if len(d) < 12 or d.OS_event.sum() < 5:
        print("  too few events — skipped"); return
    for endpt, tcol, ecol in [("OS", "OS_time", "OS_event"), ("EFS", "EFS_time", "EFS_event")]:
        dd = d.dropna(subset=[tcol, ecol])
        for score in REGS + ["composite"]:
            try:
                cph = CoxPHFitter().fit(dd[[tcol, ecol, score]], tcol, ecol)
                hr = float(np.exp(cph.params_[score])); p = float(cph.summary.loc[score, "p"])
            except Exception:
                hr, p = float("nan"), float("nan")
            hi = dd[score] > dd[score].median()
            lr = logrank_test(dd[tcol][hi], dd[tcol][~hi], dd[ecol][hi], dd[ecol][~hi])
            flag = "  *" if (p < 0.05 or lr.p_value < 0.05) else ""
            print(f"  {endpt:3s} {score:10s} Cox HR/SD={hr:5.2f} p={p:6.3f} | "
                  f"median-split log-rank p={lr.p_value:6.3f}{flag}")

analyse(df, "Pooled CBF (RUNX1-RUNX1T1 + CBFB-MYH11)")
for sb in ["RUNX1-RUNX1T1", "CBFB-MYH11"]:
    analyse(df[df.cyto == sb].reset_index(drop=True), sb)
print("\n(HR/SD = hazard ratio per 1 SD of regulon module z-score; >1 = higher score, worse outcome.)")
