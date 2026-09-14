#!/usr/bin/env python
# ====================================================================
# build_adult_scrna_umap.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_adult_scrna_umap.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Adult AML single-cell UMAP (pediatric panel-B style) — the 12 panel-5E targets on the AML
scAtlas adult-AML cells, using the atlas scVI UMAP embedding. Compartment reference + 12 gene
expression panels (atlas log-normalised X). -> assets/adult_scrna_umap.png
"""
import os
import numpy as np
import scipy.sparse as sp
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
H5 = os.path.join(ROOT, "external_data", "adult_scrna", "AML_scAtlas_748k.h5ad")
OUT = os.path.join(ROOT, "assets", "adult_scrna_umap.png")

GENES = ["CD96", "CD9", "SUCNR1", "IL1RAP", "FLT3", "CD33", "IL3RA", "CLEC12A",
         "CD7", "TNFRSF4", "ABCA7", "ITGAX"]
DISP = {"IL3RA": "CD123", "CLEC12A": "CLL-1"}
COMPART = {  # Author Cell Type -> coarse compartment
    **{c: "HSPC / progenitor" for c in ["HSPC", "GMP", "CMP", "MEP"]},
    **{c: "monocyte" for c in ["CD14+ Mono", "CD16+ Mono", "ProMono"]},
    **{c: "DC" for c in ["cDC", "pDC"]}, **{c: "T / NK" for c in ["T", "NK"]},
    **{c: "B / plasma" for c in ["B", "ProB", "Plasma"]}, "Erythroid": "erythroid"}
CORDER = ["HSPC / progenitor", "monocyte", "DC", "T / NK", "B / plasma", "erythroid"]
CCOL = {"HSPC / progenitor": "#c0392b", "monocyte": "#e67e22", "DC": "#8e44ad",
        "T / NK": "#16a085", "B / plasma": "#2980b9", "erythroid": "#7f8c8d"}

A = ad.read_h5ad(H5, backed="r")
o = A.obs
mask = ((o["Main Age Group"] == "adult") & (o["disease"] == "myeloid leukemia")).to_numpy()
idx = np.where(mask)[0]
print(f"adult AML cells for UMAP: {len(idx):,}")
emb = np.asarray(A.obsm["X_scVI_umap"])[idx]
comp = np.array([COMPART.get(str(c), "erythroid") for c in o["Author Cell Type"].values[idx]])
# malignant blasts (van Galen/Zeng leukemic hierarchy: LSC + *-like) vs NORMAL HSPC/progenitors
MALIG = {"LSC", "HSC-like", "Prog-like", "GMP-like", "cDC-like", "ProMono-like", "Mono-like"}
NORMAL_HSPC = {"HSC", "GMP", "Prog", "MPP", "CMP", "MEP"}
hct = o["HSPC Cell Type"].astype(str).values[idx]
malig = np.isin(hct, list(MALIG))
normhspc = np.isin(hct, list(NORMAL_HSPC))
print(f"  leukemic blasts: {int(malig.sum()):,} | normal HSPC: {int(normhspc.sum()):,} / {len(idx):,}")

# 12-gene expression slice (atlas log-normalised X), column-subset then row-subset
fn = A.var["feature_name"].astype(str).values
gpos = {g: int(np.where(fn == g)[0][0]) for g in GENES if (fn == g).any()}
A12 = A[:, [gpos[g] for g in GENES]].to_memory()
X = A12.X[idx]
X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)

x, y = emb[:, 0], emb[:, 1]
xl = np.percentile(x, [0.2, 99.8]); yl = np.percentile(y, [0.2, 99.8])
mpl.rcParams.update({"font.size": 9})
fig, axes = plt.subplots(3, 5, figsize=(17, 10)); axes = axes.ravel()

# (0) compartment reference (legend in the spare panel, not on the dots)
ax = axes[0]
for k in CORDER:
    m = comp == k
    ax.scatter(x[m], y[m], s=0.6, c=CCOL[k], linewidths=0, rasterized=True)
ax.set_title("compartment", fontsize=10, fontweight="bold")

# (1) leukemic blasts vs NORMAL HSPC -- where the AML cells diverge from normal progenitors
ax = axes[1]
other = ~malig & ~normhspc
ax.scatter(x[other], y[other], s=0.5, c="#e5e8ea", linewidths=0, rasterized=True)
ax.scatter(x[normhspc], y[normhspc], s=1.0, c="#2471a3", linewidths=0, rasterized=True)   # normal HSPC
ax.scatter(x[malig], y[malig], s=0.7, c="#c0392b", linewidths=0, rasterized=True)          # leukemic blasts
ax.set_title("leukemic blasts vs normal HSPC", fontsize=10, fontweight="bold")

# (2-13) gene expression
for j, g in enumerate(GENES):
    ax = axes[j + 2]
    v = X[:, j]
    vmax = max(np.percentile(v[v > 0], 98) if (v > 0).any() else 1.0, 0.5)
    ordr = np.argsort(v)
    ax.scatter(x[ordr], y[ordr], s=0.6, c=v[ordr], cmap="Reds", vmin=0, vmax=vmax,
               linewidths=0, rasterized=True)
    ax.set_title(DISP.get(g, g), fontsize=10, fontweight="bold")

ndata = len(GENES) + 2                              # 14 scatter panels (0..13)
for ax in axes[:ndata]:
    ax.set_xlim(xl); ax.set_ylim(yl); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
# spare panel (axes[14]) = legend block, off the dots
leg = axes[ndata]; leg.axis("off")
handles = [Line2D([0], [0], marker="o", linestyle="", markersize=6, markerfacecolor=CCOL[k],
           markeredgewidth=0, label=k) for k in CORDER]
handles += [Line2D([0], [0], marker="o", linestyle="", markersize=6, markerfacecolor="#c0392b",
            markeredgewidth=0, label="leukemic blast"),
            Line2D([0], [0], marker="o", linestyle="", markersize=6, markerfacecolor="#2471a3",
            markeredgewidth=0, label="normal HSPC"),
            Line2D([0], [0], marker="o", linestyle="", markersize=6, markerfacecolor="#e5e8ea",
            markeredgewidth=0, label="other normal")]
leg.legend(handles=handles, fontsize=8, loc="center", frameon=False, labelspacing=0.7,
           title="compartment / malignancy", title_fontsize=9)
for ax in axes[ndata + 1:]:
    ax.set_visible(False)
fig.suptitle("Adult AML single-cell (AML scAtlas) — 12 targets on scVI UMAP",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print("wrote", OUT)
