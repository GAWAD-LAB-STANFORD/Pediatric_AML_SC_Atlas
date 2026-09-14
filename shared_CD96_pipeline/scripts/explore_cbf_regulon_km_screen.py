#!/usr/bin/env python
# ====================================================================
# explore_cbf_regulon_km_screen.py  |  Exploratory analysis (CBF regulon survival screen)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python explore_cbf_regulon_km_screen.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""(1) KM curves for the five regulons (ERG/IKZF2/LHX6 = single-cell module z-score;
       BRCA1/E2F1 = bulk-SCENIC AUCell), median-split, OS + EFS, pooled CBF-AML.
   (2) Unbiased screen of ALL 136 bulk-SCENIC regulons for OS/EFS association in CBF-AML,
       Cox per regulon -> Benjamini-Hochberg FDR. Answers: any regulon FDR<0.05?
NO fabricated values. Reads the TARGET SCENIC loom + the two scored CSVs.
-> figures/Figure_CBF_regulon_KM.png + source_data/fig3_cbf_regulon_screen.csv
"""
import os, h5py
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests

ROOT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated"
T_LOOM = f"{ROOT}/ALSF_AML_plot/H5AD/NCI_Target_scenic_integrated-output.loom"
CLIN   = f"{ROOT}/Figure7/Figure7F_TARGET_CD96_merged.csv"
SD     = f"{ROOT}/AML_2026_Alternate/source_data"
FIG    = f"{ROOT}/AML_2026_Alternate/figures"

# ---- assemble the 5 user regulon scores (merge the two prior CSVs) + survival ----
a = pd.read_csv(f"{SD}/fig3_cbf_regulon_survival.csv")          # ERG/IKZF2/LHX6 (module z) + survival
b = pd.read_csv(f"{SD}/fig3_cbf_brca1_e2f1_survival.csv")[["sample", "BRCA1", "E2F1"]]
d5 = a.merge(b, on="sample", how="inner")
FIVE = ["ERG", "IKZF2", "LHX6", "BRCA1", "E2F1"]
print(f"CBF patients: {len(d5)} (deaths {d5.OS_event.sum()}, EFS events {d5.EFS_event.sum()})")

# ---- (1) KM figure: 2 rows (OS, EFS) x 5 regulons, median split ----
fig, axes = plt.subplots(2, 5, figsize=(20, 8))
for j, reg in enumerate(FIVE):
    for i, (ep, tc, ec) in enumerate([("OS", "OS_time", "OS_event"), ("EFS", "EFS_time", "EFS_event")]):
        ax = axes[i, j]; dd = d5.dropna(subset=[tc, ec])
        hi = dd[reg] > dd[reg].median()
        lr = logrank_test(dd[tc][hi], dd[tc][~hi], dd[ec][hi], dd[ec][~hi]).p_value
        for mask, lab, col in [(hi, "high", "#c0392b"), (~hi, "low", "#2471a3")]:
            km = KaplanMeierFitter().fit(dd[tc][mask], dd[ec][mask], label=f"{lab} (n={mask.sum()})")
            km.plot_survival_function(ax=ax, ci_show=False, color=col, linewidth=1.6)
        ax.set_title(f"{reg} — {ep}\nlog-rank p={lr:.3f}", fontsize=10,
                     fontweight="bold" if lr < 0.05 else "normal")
        ax.set_xlabel("years"); ax.set_ylim(0, 1.02); ax.legend(fontsize=7, loc="lower left")
        if j == 0: ax.set_ylabel(f"{ep} probability")
fig.suptitle("CBF-AML (RUNX1-RUNX1T1 + CBFB-MYH11, TARGET bulk, n=311) — regulon high vs low (median split)\n"
             "ERG/IKZF2/LHX6 = single-cell regulon module z-score · BRCA1/E2F1 = bulk-SCENIC AUCell",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(f"{FIG}/Figure_CBF_regulon_KM.png", dpi=160, bbox_inches="tight")
print(f"wrote {FIG}/Figure_CBF_regulon_KM.png")

# ---- (2) unbiased screen: all 136 bulk-SCENIC regulons (AUCell) ----
with h5py.File(T_LOOM, "r") as t:
    dec = lambda k: np.array([x.decode() if isinstance(x, bytes) else x for x in t["col_attrs"][k][:]])
    cid, cyto, styp = dec("CellID"), dec("Primary Cytogenetic Code"), dec("Sample.Type")
    auc = t["col_attrs"]["RegulonsAUC"]; regnames = list(auc.dtype.names)
    is_cbf = np.isin(cyto, ["RUNX1-RUNX1T1", "CBFB-MYH11"]) & np.char.startswith(styp, "Primary")
    A = pd.DataFrame({r: auc[r][:][is_cbf] for r in regnames}); A["sample"] = cid[is_cbf]
surv = d5[["sample", "OS_time", "OS_event", "EFS_time", "EFS_event"]]
A = A.merge(surv, on="sample", how="inner")
# Robust screen: median-split log-rank per regulon (the "separate into two groups" test),
# BH-FDR per endpoint over the regulons with non-zero AUCell variance.
rows = []
for ep, tc, ec in [("OS", "OS_time", "OS_event"), ("EFS", "EFS_time", "EFS_event")]:
    dd = A.dropna(subset=[tc, ec]); recs = []
    for r in regnames:
        x = dd[r].values
        if np.nanstd(x) < 1e-9:
            continue                                    # constant AUCell -> not testable
        hi = x > np.median(x)
        if hi.sum() < 5 or (~hi).sum() < 5:
            continue
        p = logrank_test(dd[tc][hi], dd[tc][~hi], dd[ec][hi], dd[ec][~hi]).p_value
        worse = "high" if dd[ec][hi].mean() > dd[ec][~hi].mean() else "low"   # group with more events
        recs.append((r.replace("_(+)", ""), ep, p, worse))
    pv = [x[2] for x in recs]; q = multipletests(pv, method="fdr_bh")[1]
    for (r, e, p, w), qq in zip(recs, q):
        rows.append(dict(regulon=r, endpoint=e, logrank_p=round(p, 4),
                         q_BH=round(float(qq), 3), worse_group=w))
res = pd.DataFrame(rows).sort_values("logrank_p")
res.to_csv(f"{SD}/fig3_cbf_regulon_screen.csv", index=False)
print(f"\n=== Unbiased screen: 136 bulk regulons x (OS,EFS) in CBF-AML, BH-FDR ===")
print("Top 12 by raw p:"); print(res.head(12).to_string(index=False))
sig = res[res.q_BH < 0.05]
print(f"\nRegulons with FDR q<0.05: {len(sig)}")
print(sig.to_string(index=False) if len(sig) else "  NONE — no regulon separates CBF-AML survival after multiple-testing correction.")
print(f"\nmin BH-q (OS) = {res[res.endpoint=='OS'].q_BH.min():.3f} ; min BH-q (EFS) = {res[res.endpoint=='EFS'].q_BH.min():.3f}")
