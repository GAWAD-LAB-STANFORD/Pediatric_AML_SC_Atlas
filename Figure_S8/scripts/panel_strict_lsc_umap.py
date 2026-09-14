#!/usr/bin/env python
"""Supplementary panel B: the strict true-LSC candidates on the atlas UMAP (leukemic-
localized, not the normal-HSC cluster). Reads strict_LSC_barcodes.csv (compute_strict_lsc.py)
and the atlas Cell Type / UMAP. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
OUT = os.path.join(C.BASE, "04_Main_Figures/Figure_2_panels")
strict = set(pd.read_csv(os.path.join(C.DATA, "strict_LSC_barcodes.csv"))["barcode"].astype(str))
f = h5py.File(C.H5, "r")
names = np.array([x.decode() if isinstance(x, bytes) else x for x in f["obs"]["_index"][:]])
umap = f["obsm/X_umap"][:]
nd = f["obs"]["Cell Type"]; cc = [x.decode() if isinstance(x, bytes) else x for x in nd["categories"][:]]
ctype = np.array([cc[i] if i >= 0 else "NA" for i in nd["codes"][:]], dtype=object); f.close()
c = np.array(["background"] * len(names), dtype=object)
c[ctype == "0_HSPC"] = "normal HSPC"
c[np.array([b in strict for b in names])] = "strict true-LSC"
fig, ax = plt.subplots(figsize=(6.4, 5.8))
for k, s, col, al in [("background", 2, "#ececec", 0.4), ("normal HSPC", 7, "#94D2BD", 0.9)]:
    m = c == k; ax.scatter(umap[m, 0], umap[m, 1], s=s, c=col, linewidths=0, rasterized=True, label=(k if k != "background" else None), alpha=al)
m = c == "strict true-LSC"
ax.scatter(umap[m, 0], umap[m, 1], s=90, c="#9B2226", edgecolors="black", linewidths=1.1, label=f"strict true-LSC (n={int(m.sum())})", zorder=5)
ax.set_xticks([]); ax.set_yticks([]); ax.set_frame_on(False)
ax.legend(loc="lower left", frameon=False, fontsize=10, markerscale=1.1)
ax.set_title("Strict LSC (CD34+CD38- , HLF/AVP+ quiescent, non-cycling, TIM3/CD96/CLEC12A+):\n%d cells, 2 patients" % int(m.sum()), fontsize=9)
plt.tight_layout()
for e in ("png", "pdf"): fig.savefig(f"{OUT}/UMAP_strict_LSC.{e}", dpi=200, bbox_inches="tight")
print(f"wrote UMAP_strict_LSC (n={int(m.sum())})")
