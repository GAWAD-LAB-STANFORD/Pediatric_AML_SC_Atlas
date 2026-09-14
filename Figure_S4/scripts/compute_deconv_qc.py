#!/usr/bin/env python
"""S4 (bulk-deconvolution validation) DATA for the corrected 49-state scheme. From the deployed
CIBERSORTx deconvolution of TARGET bulk (LS_survival.csv = Jobs 24-27, 49 leukemic states + 8
normal compartments): (1) per-sample goodness-of-fit (CIBERSORTx Pearson correlation, RMSE,
P-value) across all deconvolved samples; (2) per-49-LS detection reliability (fraction of samples
with fraction >1% and >5%, and mean fraction). Replaces the legacy AML_1..28 pseudobulk figure and
the stale 30-state detection table. No fabricated values."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"); import config as C
D2 = C.DATA.replace("Figure_3", "Figure_2")
d = pd.read_csv(os.path.join(C.DATA, "LS_survival.csv"))
n = len(d)
# per-sample fit QC
qc = d[["sample", "P-value", "Correlation", "RMSE"]].copy()
for c in ["P-value", "Correlation", "RMSE"]: qc[c] = pd.to_numeric(qc[c], errors="coerce")
qc.to_csv(os.path.join(D2, "deconv_fit_qc.csv"), index=False)
# per-49-LS detection across all deconvolved TARGET samples
ls = [c for c in d.columns if c.startswith("LS_")]
F = d[ls].apply(pd.to_numeric, errors="coerce")
det = pd.DataFrame({"LS": ls,
                    "det_gt1pct": (F > 0.01).mean().values,
                    "det_gt5pct": (F > 0.05).mean().values,
                    "mean_frac": F.mean().values}).round(4)
det.to_csv(os.path.join(D2, "LS_detection_49.csv"), index=False)
print(f"n deconvolved TARGET samples = {n}")
print(f"fit: median Correlation={qc.Correlation.median():.3f}, median RMSE={qc.RMSE.median():.3f}, "
      f"frac P<0.05={(qc['P-value']<0.05).mean():.3f}")
print(f"detection: median det>1% = {det.det_gt1pct.median():.3f}; "
      f"states det>1%<0.20 = {det.LS[det.det_gt1pct<0.20].tolist()}")
print("wrote deconv_fit_qc.csv, LS_detection_49.csv")
