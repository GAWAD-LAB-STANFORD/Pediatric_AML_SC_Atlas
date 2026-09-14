#!/usr/bin/env python
# ====================================================================
# rebuild_figure1_umaps.py  |  CD96 figure pipeline component
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python rebuild_figure1_umaps.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Rebuild Figure 1 UMAP panels from the source h5ad.

Reads only obsm/X_umap + the relevant obs columns + stored colour palettes
from ALSF_AML_Combo_3500_with_PAC_new.h5ad (no expression matrix loaded),
builds a lightweight AnnData, and regenerates the UMAP panels with scanpy
to match the committed figures.
"""
import os
import h5py
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
H5AD = os.path.join(HERE, "..", "ALSF_AML_plot", "H5AD",
                    "ALSF_AML_Combo_3500_with_PAC_new.h5ad")
OUTDIR = os.path.join(HERE, "rebuilt")
os.makedirs(OUTDIR, exist_ok=True)

OBS_COLS = ["Cell Type", "lineage", "Prognosis-Associated Clusters",
            "FAB", "Cytogenetic"]

# The h5ad's stored `lineage_colors` is the 28-colour palette, but the
# published Figure 1 lineage panel uses scanpy's 20-colour palette (the same
# one the committed FAB/Cytogenetic/PAC panels use). Override to match.
DEFAULT_20 = ["#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc",
              "#8c564b", "#e377c2", "#b5bd61", "#17becf", "#aec7e8",
              "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94",
              "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31"]
PALETTE_OVERRIDE = {"lineage": DEFAULT_20}


def _decode(arr):
    return [x.decode() if isinstance(x, bytes) else x for x in arr]


def read_categorical(f, col):
    g = f["obs"][col]
    cats = _decode(g["categories"][:])
    codes = g["codes"][:].astype(int)
    series = pd.Categorical.from_codes(codes, categories=cats)
    colors = None
    ck = col + "_colors"
    if ck in f["uns"]:
        stored = _decode(f["uns"][ck][:])
        if len(stored) >= len(cats):
            colors = stored[:len(cats)]   # first N align with current categories
    return series, colors


def main():
    print(f"Reading {os.path.basename(H5AD)} (metadata only) ...")
    with h5py.File(H5AD, "r") as f:
        umap = f["obsm/X_umap"][:]
        obs = {}
        palettes = {}
        for col in OBS_COLS:
            series, colors = read_categorical(f, col)
            obs[col] = series
            if col in PALETTE_OVERRIDE:
                colors = PALETTE_OVERRIDE[col][:len(series.categories)]
            if colors is not None:
                palettes[col] = colors
    n = umap.shape[0]
    print(f"  {n} cells")

    adata = ad.AnnData(
        X=np.zeros((n, 1), dtype=np.float32),
        obs=pd.DataFrame(obs),
        var=pd.DataFrame(index=["_dummy"]),
    )
    adata.obsm["X_umap"] = umap
    for col, colors in palettes.items():
        adata.uns[col + "_colors"] = list(colors)

    sc.settings.figdir = OUTDIR
    sc.set_figure_params(dpi=150, dpi_save=300, frameon=False)

    panels = [
        ("Cell Type",                     "Cell Type",                     "umap_ALSF_AML_Cell_Type.pdf"),
        ("lineage",                       "Lineage",                       "umap_ALSF_AML_Lineage.pdf"),
        ("Prognosis-Associated Clusters", "Prognosis-Associated Clusters", "umap_PAC_clusters.pdf"),
    ]
    for col, title, fname in panels:
        sc.pl.umap(adata, color=col, title=title, frameon=False, show=False)
        out = os.path.join(OUTDIR, fname)
        plt.savefig(out, bbox_inches="tight")
        plt.close()
        print(f"  wrote {fname}")

    # FAB + Cytogenetic as a combined two-panel PNG
    sc.pl.umap(adata, color=["FAB", "Cytogenetic"], frameon=False, show=False)
    out = os.path.join(OUTDIR, "umap_Combo_FAB_CYtogenetic.png")
    plt.savefig(out, bbox_inches="tight", dpi=300)
    plt.close()
    print("  wrote umap_Combo_FAB_CYtogenetic.png")
    print(f"Done -> {OUTDIR}")


if __name__ == "__main__":
    main()
