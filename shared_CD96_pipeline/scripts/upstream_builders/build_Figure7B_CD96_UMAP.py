#!/usr/bin/env python
# ====================================================================
# build_Figure7B_CD96_UMAP.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_Figure7B_CD96_UMAP.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Figure 7B -- scRNA atlas: Cell Type compartments + CD96 expression overlay.

Two side-by-side UMAPs on the same coordinates:
  left   compartment-grouped Cell Type   (8 functional groups)
  right  CD96 expression                  (log-normalized, Reds colormap)

The 29 Cell Type categories in the h5ad include the AML-prefix labels that
identify normal cells captured *from* AML samples (these are transcriptionally
normal T / NK / B / monocyte / erythroid cells -- see project memory).  For
visual clarity, all real leukemic clusters (AML, AML-MKI67, AML-PCNA) are
collapsed into one "AML (leukemic)" group; the AML-prefixed normal cells
are merged with their normal counterparts; HSPC and Myeloid_Pro are kept
separate because they anchor the safety story.
"""
import os
import sys
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, "supplements")
sys.path.insert(0, SUPP)
from _fig7_style import PALETTE, RC, save_fig, SEQ_CMAP, CD96_CMAP

H5AD = os.path.join(HERE, "..", "ALSF_AML_plot", "H5AD",
                    "ALSF_AML_Combo_3500_with_PAC_new.h5ad")

# 29 Cell Type categories -> 8 functional compartments
CELLTYPE_GROUPING = {
    # leukemic blasts (real AML)
    "AML":                "AML (leukemic)",
    "AML-MKI67":          "AML (leukemic)",
    "AML-PCNA":           "AML (leukemic)",
    # safety anchors
    "0_HSPC":             "HSPC",
    "Myeloid_Pro":        "Myeloid_Pro",
    # T cells  (AML-prefixed labels are normal T cells from AML samples)
    "CTL":                "T cells",
    "Naïve_CD4T":         "T cells",
    "Naïve_CD8T":         "T cells",
    "Activated_CD4T":     "T cells",
    "AML-CTL":            "T cells",
    "AML-CD4T":           "T cells",
    "AML-Naïve_CD8T":     "T cells",
    # NK
    "NK":                 "NK",
    "AML-NK":             "NK",
    # B-lineage (pro-B, pre-B, mature B, plasma)
    "CD34+ProB":          "B-lineage",
    "ProB":               "B-lineage",
    "PreB":               "B-lineage",
    "CD20+B":             "B-lineage",
    "PlasmaB":            "B-lineage",
    "AML-B":              "B-lineage",
    # mature myeloid (monocytes, macrophages, DCs)
    "CD14_Monocyte":      "Mature myeloid",
    "CD16_Monocyte":      "Mature myeloid",
    "Macrophage":         "Mature myeloid",
    "mDC":                "Mature myeloid",
    "pDC":                "Mature myeloid",
    "AML-CD14":           "Mature myeloid",
    "AML-CD1C":           "Mature myeloid",
    # erythroid
    "Erythrocytes":       "Erythroid",
    "AML-Ery":            "Erythroid",
}

# canonical compartment ordering + colors
GROUP_ORDER = ["AML (leukemic)", "HSPC", "Myeloid_Pro",
                "Mature myeloid", "B-lineage", "T cells", "NK", "Erythroid"]
GROUP_COLORS = {
    "AML (leukemic)":   PALETTE["CD96"],     # crimson  -- emphasize the target compartment
    "HSPC":             "#27ae60",            # green    -- safety anchor (HSPC-sparing)
    "Myeloid_Pro":      "#16a085",            # teal     -- safety anchor
    "Mature myeloid":   "#bdc3c7",            # light grey
    "B-lineage":        "#2980b9",            # blue
    "T cells":          PALETTE["T_cell"],    # purple
    "NK":               "#9b59b6",            # lighter purple
    "Erythroid":        "#7f8c8d",            # grey
}


def _decode(arr):
    return [x.decode() if isinstance(x, bytes) else x for x in arr]


def main():
    print("  loading h5ad ...")
    with h5py.File(H5AD, "r") as f:
        # cell type
        ctype_cats = _decode(f["obs"]["Cell Type"]["categories"][:])
        ctype = np.array(ctype_cats)[f["obs"]["Cell Type"]["codes"][:].astype(int)]
        # UMAP
        umap = f["obsm/X_umap"][:]
        # raw CD96 counts (sparse CSR) -> log-normalized
        rx = f["raw"]["X"]
        X = sp.csr_matrix((rx["data"][:], rx["indices"][:], rx["indptr"][:]),
                          shape=tuple(rx.attrs["shape"])).tocsc()
        rk = f["raw"]["var"].attrs.get("_index", b"_index")
        rk = rk.decode() if isinstance(rk, bytes) else rk
        symbols = _decode(f["raw"]["var"][rk][:])
        cd96_idx = symbols.index("CD96")
        cd96_raw = np.asarray(X[:, cd96_idx].todense(), dtype=np.float32).ravel()
        # library size for normalization
        total = np.asarray(X.sum(axis=1)).ravel().astype(np.float32)

    n = len(ctype)
    print(f"    cells: {n:,}")
    print(f"    UMAP: {umap.shape}")

    # log-normalize CD96 the standard scanpy way
    cd96_norm = np.log1p(cd96_raw / np.maximum(total, 1) * 10_000).astype(np.float32)
    print(f"    CD96 log-norm  min={cd96_norm.min():.2f}  "
          f"max={cd96_norm.max():.2f}  "
          f"%>0={(cd96_norm > 0).mean()*100:.1f}")

    # group cell types
    unmapped = sorted(set(ctype) - set(CELLTYPE_GROUPING))
    if unmapped:
        print(f"  WARNING unmapped Cell Type categories (assigned to 'Other'): {unmapped}")
    group = np.array([CELLTYPE_GROUPING.get(c, "Other") for c in ctype])

    # build AnnData (CD96 in X so sc.pl.umap can paint by it)
    adata = ad.AnnData(
        X=cd96_norm.reshape(-1, 1),
        obs=pd.DataFrame({
            "compartment": pd.Categorical(group, categories=GROUP_ORDER),
        }),
        var=pd.DataFrame(index=["CD96"]),
    )
    adata.obsm["X_umap"] = umap
    adata.uns["compartment_colors"] = [GROUP_COLORS[g] for g in GROUP_ORDER]

    # --- render ---
    with mpl.rc_context(RC):
        sc.set_figure_params(dpi=150, dpi_save=300, frameon=False)
        fig, axes = plt.subplots(1, 2, figsize=(15, 6.6),
                                  gridspec_kw={"wspace": 0.32,
                                                "top": 0.90,
                                                "bottom": 0.04,
                                                "left": 0.02,
                                                "right": 0.98})

        # left: compartments  -- suppress scanpy's own title; we set ours below
        sc.pl.umap(adata, color="compartment", ax=axes[0],
                    show=False, frameon=False, size=8,
                    legend_loc="right margin", title="",
                    legend_fontsize=9)
        axes[0].set_title("cell-type compartment",
                          loc="left", fontsize=10.5, fontweight="bold",
                          pad=8)

        # right: CD96 expression  -- CD96 RED colormap (its identity colour), so
        # the lead target reads red throughout the figure; vmax clipped so the
        # mid-expression bulk of the AML compartment is visible.
        sc.pl.umap(adata, color="CD96", ax=axes[1],
                    show=False, frameon=False, size=8,
                    cmap=CD96_CMAP, vmin=0, vmax=2.5, title="")
        axes[1].set_title("CD96 expression", loc="left", fontsize=10.5,
                          fontweight="bold", pad=8)
        # (figure title/commentary in caption -- added by authors)
        save_fig(fig, HERE, "Figure7B_CD96_UMAP")
    print("\nDone -> Figure7B_CD96_UMAP.{pdf,png}")


if __name__ == "__main__":
    main()
