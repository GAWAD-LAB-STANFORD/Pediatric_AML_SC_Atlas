#!/usr/bin/env python
"""Corrected 'Relapsed' covariate overlay for S2 (AML atlas covariates). The h5ad obs['Relapsed']
is unreliable (memory relapse-annotation-source: AML3121 mislabeled, AML4363 unknown). This rebuilds
the atlas UMAP coloured by the AUTHORITATIVE per-sample relapse status reconciled from the clinical
xlsx (MRD_BCB-AML-list, 'AML list'). Corrected relapsed set = 8 samples; AML4271/AML4363 unknown.
Replaces the stale sub-panel. No fabricated values."""
import os, numpy as np, pandas as pd, scanpy as sc, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
H5 = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
OUT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/05_Supplementary_Figures"
RELAPSED = {"AML948","AML3082","AML3371","AML3492","AML4068","AML4239","AML4304","AML3121"}  # authoritative
UNKNOWN  = {"AML4271","AML4363"}                                                            # clinically unavailable/ambiguous
a = sc.read_h5ad(H5, backed="r")
um = np.asarray(a.obsm["X_umap"]); sid = a.obs["SampleID"].astype(str).values
def cls(s):
    if s.startswith("0_HealthyBM"): return "healthy BM"
    if s in RELAPSED: return "relapsed"
    if s in UNKNOWN:  return "unknown"
    return "not relapsed"
cat = np.array([cls(s) for s in sid])
col = {"not relapsed":"#c7ccd1","relapsed":"#d1495b","unknown":"#8d99ae","healthy BM":"#2a9d8f"}
order = ["not relapsed","unknown","healthy BM","relapsed"]   # draw relapsed last (on top)
fig, ax = plt.subplots(figsize=(6,5.4))
for c in order:
    m = cat == c
    ax.scatter(um[m,0], um[m,1], s=2, c=col[c], label=f"{c} (n={m.sum()})", linewidths=0, rasterized=True)
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
ax.set_title("Relapse status (corrected from clinical record)", fontsize=11)
lg = ax.legend(loc="lower left", fontsize=7, markerscale=3, frameon=False)
for e in ("png","pdf"): fig.savefig(os.path.join(OUT, f"Figure_S2_panelB_relapsed_corrected.{e}"), dpi=200, bbox_inches="tight")
print("relapsed cells:", int((cat=='relapsed').sum()), "| samples relapsed:", sorted(RELAPSED))
print("wrote Figure_S2_panelB_relapsed_corrected")
