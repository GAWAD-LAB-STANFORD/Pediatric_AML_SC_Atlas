#!/usr/bin/env python
"""Figure 3B (LS19): LS19 (poor) with normal pDC (mature reference) and normal HSPC (progenitor
reference) on the atlas UMAP. Shows LS19 sits with the progenitor/HSPC compartment, distinct from
mature pDC. Data = h5ad UMAP + REFA/Cell_Type. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"); import config as C
OUT = os.path.join(C.BASE, "04_Main_Figures/Figure_3_panels")
f = h5py.File(C.H5, "r"); dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
obs = dec(f["obs"]["_index"][:]); um = f["obsm"]["X_umap"][:]
ct = pd.Series(dec(f["obs"]["Cell_Type"]["categories"][:])[f["obs"]["Cell_Type"]["codes"][:]], index=obs); f.close()
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
def grp(b):
    if bc2ls.get(b) == "LS_19": return "LS19 (poor)"
    if ct.get(b) == "pDC": return "normal pDC"
    if ct.get(b) == "0_HSPC": return "normal HSPC"
    return "background"
g = np.array([grp(b) for b in obs])
# colours drawn from the shared manuscript SEQ palette (the dot-plot scheme): teal progenitor,
# amber mature pDC, deep-red LS19
col = {"background": "#e9e9e9", "normal HSPC": "#0A9396", "normal pDC": "#EE9B00", "LS19 (poor)": "#AE2012"}
order = ["background", "normal HSPC", "normal pDC", "LS19 (poor)"]
fig, ax = plt.subplots(figsize=(6, 5.4))
for c in order:
    m = g == c
    ax.scatter(um[m, 0], um[m, 1], s=(2 if c == "background" else 5), c=col[c], linewidths=0,
               rasterized=True, label=f"{c} (n={m.sum()})" if c != "background" else None)
ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
ax.legend(loc="lower left", fontsize=8, markerscale=3, frameon=False)
for e in ("png", "pdf"): fig.savefig(os.path.join(OUT, f"Figure_3_LS19_umap.{e}"), dpi=200, bbox_inches="tight")
print("wrote Figure_3_LS19_umap")
