#!/usr/bin/env python
"""Figure 3B: LS10 (poor) vs the non-prognostic proliferation comparators (LS39, LS3) and
gated normal HSC on the atlas UMAP. LS10 sits nearest HSC; LS3 drifts toward myeloid.
Cells/labels from reference_assignments + the gated-HSC barcode list. No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import config as C
D2 = C.DATA.replace("Figure_3", "Figure_2")
OUT = os.path.join(C.BASE, "04_Main_Figures/Figure_2_panels")
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}; ra["LS"] = ra["res4.50"].map(raw2ls)
lab = {"LS_10": "LS10 (poor)", "LS_39": "LS39 (ns)", "LS_3": "LS3 (ns)"}
bc2 = {b: lab[s] for b, s in zip(ra.barcode, ra.LS) if s in lab}
gated = set(pd.read_csv(os.path.join(D2, "normal_HSC_gated_barcodes.csv"))["barcode"].astype(str))
f = h5py.File(C.H5, "r")
names = [x.decode() if isinstance(x, bytes) else x for x in f["obs"]["_index"][:]]
umap = f["obsm/X_umap"][:]; f.close()
cat = np.array(["background"] * len(names), dtype=object)
for i, b in enumerate(names):
    if b in gated: cat[i] = "normal HSC"
    elif b in bc2: cat[i] = bc2[b]
order = ["background", "LS39 (ns)", "LS3 (ns)", "LS10 (poor)", "normal HSC"]
col = {"background": "#eeeeee", "LS39 (ns)": "#0072B2", "LS3 (ns)": "#009E73", "LS10 (poor)": "#d1495b", "normal HSC": "#E69F00"}
fig, ax = plt.subplots(figsize=(6.4, 5.8))
for k in order:
    m = cat == k
    ax.scatter(umap[m, 0], umap[m, 1], s=(1.2 if k == "background" else (9 if k == "normal HSC" else (7 if k == "LS10 (poor)" else 3.5))),
               c=col[k], linewidths=0, rasterized=True, label=k, alpha=(0.3 if k == "background" else 0.9))
ax.set_xticks([]); ax.set_yticks([]); ax.set_frame_on(False)
ax.legend(loc="lower left", frameon=False, markerscale=3, fontsize=9, handletextpad=0.2)
plt.tight_layout()
for e in ("png", "pdf"): fig.savefig(f"{OUT}/UMAP_LS10_vs_prolifNS_HSC.{e}", dpi=200, bbox_inches="tight")
print("wrote UMAP_LS10_vs_prolifNS_HSC")
