#!/usr/bin/env python
# ====================================================================
# build_fig5B_umaps.py  |  Figure 5 (builder)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_fig5B_umaps.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""ALTERNATE Figure 5 panel B — UMAPs for the 4 LEAD targets (CD96 + CD9, SUCNR1, IL1RAP),
plus a cell-type-compartment reference, from the bundled count matrix's X_umap embedding.
Embedded as an image (the agreed scanpy/UMAP exception). -> assets/fig6B_leads_umap.png
"""
import os
import numpy as np
import scipy.sparse as sp
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

# manuscript heatmap0 sequential colormap (matches Figure 3 scale_fill_seq); low->dark teal, high->warm red
SEQ_CMAP = LinearSegmentedColormap.from_list("heatmap0",
    ["#001219", "#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#AE2012", "#9B2226"])

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
CM = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
OUT = os.path.join(ROOT, "assets", "fig6B_leads_umap.png")
LEADS = ["CD96", "CD9", "SUCNR1", "IL1RAP"]   # AMN dropped (4 leads)
FLAG = {}   # no dagger; CD9 shown plain (caveat lives in the safety panel)
LEAD_COL = {"CD96": "#c0392b", "CD9": "#0072B2", "SUCNR1": "#009E73", "AMN": "#E69F00", "IL1RAP": "#56B4E9"}

COMPART = {  # cell type -> coarse compartment (for the reference UMAP)
    **{c: "leukemic" for c in ["AML", "AML-MKI67", "AML-PCNA"]},
    **{c: "HSPC / Myeloid" for c in ["0_HSPC", "Myeloid_Pro"]},
    **{c: "monocyte / DC" for c in ["CD14_Monocyte", "CD16_Monocyte", "mDC", "pDC", "Macrophage", "AML-CD14", "AML-CD1C"]},
    **{c: "T / NK" for c in ["Naïve_CD4T", "Naïve_CD8T", "CTL", "Activated_CD4T", "NK", "AML-CD4T", "AML-CTL", "AML-NK", "AML-Naïve_CD8T"]},
    **{c: "B / plasma" for c in ["CD20+B", "ProB", "PreB", "CD34+ProB", "PlasmaB", "AML-B"]},
    **{c: "erythroid" for c in ["Erythrocytes", "AML-Ery"]},
}
COMP_ORDER = ["leukemic", "HSPC / Myeloid", "monocyte / DC", "T / NK", "B / plasma", "erythroid"]
COMP_COL = {"leukemic": "#c0392b", "HSPC / Myeloid": "#e67e22", "monocyte / DC": "#8e44ad",
            "T / NK": "#16a085", "B / plasma": "#2980b9", "erythroid": "#7f8c8d"}

A = ad.read_h5ad(CM)
um = A.obsm["X_umap"]
ct = A.obs["Cell Type"].astype(str).values
comp = np.array([COMPART.get(c, "erythroid") for c in ct])
# expression = the atlas log-normalized 'counts' layer (the same layer used for the panel-D dot plots)
Xl = A[:, LEADS].layers["counts"]
expr = np.asarray(Xl.todense()) if sp.issparse(Xl) else np.asarray(Xl)

x, y = um[:, 0], um[:, 1]
xlim = np.percentile(x, [0.2, 99.8]); ylim = np.percentile(y, [0.2, 99.8])
order = np.argsort(np.random.RandomState(0).rand(len(x)))  # stable shuffle for overplot

mpl.rcParams.update({"font.size": 9})
fig, axes = plt.subplots(2, 3, figsize=(11, 7.2))
axes = axes.ravel()

# (0) cell-type compartment reference
ax = axes[0]
for k in COMP_ORDER:
    m = comp == k
    ax.scatter(x[m], y[m], s=1.2, c=COMP_COL[k], linewidths=0, rasterized=True)
ax.set_title("cell compartment", fontsize=10, fontweight="bold")
ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="", markersize=4, markerfacecolor=COMP_COL[k],
                          markeredgewidth=0, label=k) for k in COMP_ORDER],
          fontsize=6, loc="lower left", frameon=False, handletextpad=0.1, labelspacing=0.25)

# (1-5) expression UMAPs
for j, g in enumerate(LEADS):
    ax = axes[j + 1]
    v = expr[:, j]
    vmax = max(np.percentile(v[v > 0], 98) if (v > 0).any() else 1.0, 0.5)
    o = np.argsort(v)  # plot high-expressers on top
    sc = ax.scatter(x[o], y[o], s=1.2, c=v[o], cmap=SEQ_CMAP, vmin=0, vmax=vmax,
                    linewidths=0, rasterized=True)
    ttl = g + FLAG.get(g, "")
    ax.set_title(ttl, fontsize=11, fontweight="bold", color=LEAD_COL[g])
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.ax.tick_params(labelsize=6); cb.outline.set_visible(False)

for ax in axes[len(LEADS) + 1:]:        # hide unused cells (4 genes + 1 ref = 5 of 6)
    ax.set_visible(False)
for ax in axes[:len(LEADS) + 1]:
    ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
# no figure suptitle; commentary -> caption (per-UMAP gene labels + 'cell compartment' kept)
fig.tight_layout(rect=[0, 0.0, 1, 1.0])
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
