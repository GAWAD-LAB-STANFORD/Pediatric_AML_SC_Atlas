#!/usr/bin/env python
# ====================================================================
# build_Figure7D_cyto_PAC_heatmap.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_Figure7D_cyto_PAC_heatmap.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Figure 6D -- AML surface-target coverage across cytogenetic subtypes and
PAC clusters (scRNA).  Nine targets x subtypes/clusters; colour = % cells with
>=1 raw read.

Two stacked heatmaps, shared colourmap:
  TOP    -- cytogenetic subtypes (leukemic cells only), columns sorted by CD96
            coverage.  Shows CD96 covers most subtypes but is cold in
            KMT2Ar(MLLr) / NUP98-NSD1 / MYB-GATA1 / del7q -- and how the other
            targets fill (or share) those gaps.
  BOTTOM -- the 9 PAC clusters (all cells), biological order.  HSPC +
            Myeloid_Pro stay clean for CD96 (safety); the comparators do not.

Pure data from the 28-patient atlas; "% positive" = raw count > 0.
"""
import os
import sys
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, "supplements")
sys.path.insert(0, SUPP)
from _fig7_style import PALETTE, RC, save_fig, SEQ_CMAP, CD96_CMAP

H5AD = os.path.join(HERE, "..", "ALSF_AML_plot", "H5AD",
                    "ALSF_AML_Combo_3500_with_PAC_new.h5ad")

GENES = ["CD96", "CD33", "IL3RA", "CLEC12A", "FLT3", "MSLN",
         "CD70", "ADGRG1", "IL1RAP"]
GLABEL = {
    "CD96": "CD96", "CD33": "CD33", "IL3RA": "CD123 (IL3RA)",
    "CLEC12A": "CLL-1 (CLEC12A)", "FLT3": "FLT3", "MSLN": "MSLN",
    "CD70": "CD70", "ADGRG1": "GPR56 (ADGRG1)", "IL1RAP": "IL1RAP",
}

LEUK_SUBS = [f"AML_{i}" for i in range(1, 28)]
PAC_ORDER = ["0_HSPC", "Myeloid_Pro", "FPAC_1", "FPAC_2",
             "PPAC_1", "PPAC_2", "PPAC_3", "PPAC_4", "PPAC_5"]
SAFETY_PACS = {"0_HSPC", "Myeloid_Pro"}
ADVERSE_PPACS = {"PPAC_1", "PPAC_2", "PPAC_3", "PPAC_4", "PPAC_5"}
FAVORABLE_PACS = {"FPAC_1", "FPAC_2"}
MIN_CELLS_PER_CYTO = 500
VMAX = 70.0


def _decode(arr):
    return [x.decode() if isinstance(x, bytes) else x for x in arr]


def _cat(f, key):
    return np.array(pd.Categorical.from_codes(
        f["obs"][key]["codes"][:].astype(int),
        categories=_decode(f["obs"][key]["categories"][:])).astype(object))


def main():
    print("  loading h5ad ...")
    with h5py.File(H5AD, "r") as f:
        cyto = _cat(f, "Cytogenetic")
        pac = _cat(f, "Prognosis-Associated Clusters")
        sub = _cat(f, "AML Sub-Clusters")
        rx = f["raw"]["X"]
        X = sp.csr_matrix((rx["data"][:], rx["indices"][:], rx["indptr"][:]),
                          shape=tuple(rx.attrs["shape"])).tocsc()
        rk = f["raw"]["var"].attrs.get("_index", b"_index")
        rk = rk.decode() if isinstance(rk, bytes) else rk
        symbols = _decode(f["raw"]["var"][rk][:])
        pos = {}
        for g in GENES:
            if g not in symbols:
                raise KeyError(f"{g} not in h5ad raw var")
            pos[g] = (np.asarray(X[:, symbols.index(g)].todense(),
                                 np.float32).ravel() > 0)

    leuk = np.isin(sub, LEUK_SUBS)

    # cytogenetic columns (leukemic cells, n>=MIN), order by CD96 desc
    cyto_cols, cyto_n = [], {}
    for c in sorted(set(cyto[leuk])):
        if "HealthyBM" in str(c):
            continue
        m = leuk & (cyto == c)
        if int(m.sum()) >= MIN_CELLS_PER_CYTO:
            cyto_cols.append(c); cyto_n[c] = int(m.sum())
    cyto_cols.sort(key=lambda c: pos["CD96"][leuk & (cyto == c)].mean(),
                   reverse=True)
    Mc = np.array([[100.0 * pos[g][leuk & (cyto == c)].mean()
                    for c in cyto_cols] for g in GENES])

    # PAC columns (all cells), biological order
    pac_cols = [p for p in PAC_ORDER if (pac == p).any()]
    pac_n = {p: int((pac == p).sum()) for p in pac_cols}
    Mp = np.array([[100.0 * pos[g][pac == p].mean() for p in pac_cols]
                   for g in GENES])

    pd.DataFrame(Mc, index=GENES, columns=cyto_cols).round(3).to_csv(
        os.path.join(HERE, "Figure7D_target_by_cyto.csv"))
    pd.DataFrame(Mp, index=GENES, columns=pac_cols).round(3).to_csv(
        os.path.join(HERE, "Figure7D_target_by_PAC.csv"))
    print("\n  --- CD96 row, cytogenetic (sorted) ---")
    print(pd.Series(Mc[0], index=cyto_cols).round(1).to_string())

    cmap = SEQ_CMAP
    norm = mpl.colors.Normalize(vmin=0, vmax=VMAX)

    def s_of(v):              # smaller scale: denser grid, avoid dot overlap
        return v * 4.5 + 1.5

    def draw(ax, M, cols, ncount, kind):
        ng2, ncols = len(GENES), len(cols)
        if kind == "pac":
            for j, name in enumerate(cols):
                if name in SAFETY_PACS:
                    ax.axvspan(j - 0.5, j + 0.5, color="#eaf5ee", zorder=0)
        ax.axhspan(-0.5, 0.5, color="#fbecea", zorder=0)   # CD96 row
        ax.set_axisbelow(True)
        ax.grid(color="#e9eef4", lw=0.6, zorder=1)
        # other targets -> navy magnitude colormap;  CD96 row (i==0) -> red
        xs, ys, ss, cc = [], [], [], []
        rx, ry, rs, rc = [], [], [], []
        for i in range(ng2):
            for j in range(ncols):
                if i == 0:
                    rx.append(j); ry.append(i); rs.append(s_of(M[i, j])); rc.append(M[i, j])
                else:
                    xs.append(j); ys.append(i); ss.append(s_of(M[i, j])); cc.append(M[i, j])
        ax.scatter(xs, ys, s=ss, c=cc, cmap=cmap, norm=norm,
                   edgecolor="#3b5870", linewidth=0.35, zorder=3)
        ax.scatter(rx, ry, s=rs, c=rc, cmap=CD96_CMAP, norm=norm,
                   edgecolor="#7d2018", linewidth=0.45, zorder=4)
        ax.set_xlim(-0.6, ncols - 0.4); ax.set_ylim(-0.6, ng2 - 0.4)
        ax.invert_yaxis()
        ax.set_yticks(range(ng2))
        yt = ax.set_yticklabels([GLABEL[g] for g in GENES], fontsize=8.5)
        yt[0].set_color(PALETTE["CD96"]); yt[0].set_fontweight("bold")
        ax.set_xticks(range(ncols))
        xt = ax.set_xticklabels([f"{c}\n(n={ncount[c]:,})" for c in cols],
                                rotation=35, ha="right", fontsize=7.8)
        if kind == "pac":
            for lab, name in zip(xt, cols):
                if name in SAFETY_PACS:
                    lab.set_color(PALETTE["favorable"]); lab.set_fontweight("bold")
                elif name in ADVERSE_PPACS:
                    lab.set_color(PALETTE["adverse"]); lab.set_fontweight("bold")
                elif name in FAVORABLE_PACS:
                    lab.set_color("#3498db")
        else:
            for lab, name in zip(xt, cols):
                if M[0, cols.index(name)] < 10:   # CD96-cold subtype
                    lab.set_color("#7f8c8d"); lab.set_fontweight("bold")
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)

    with mpl.rc_context(RC):
        fig = plt.figure(figsize=(14.0, 8.8))
        gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 0.12],
                              height_ratios=[1.0, 1.0], wspace=0.03, hspace=0.6)
        ax_cy = fig.add_subplot(gs[0, 0])
        ax_pa = fig.add_subplot(gs[1, 0])
        ax_leg = fig.add_subplot(gs[:, 1]); ax_leg.axis("off")
        draw(ax_cy, Mc, cyto_cols, cyto_n, "cyto")
        draw(ax_pa, Mp, pac_cols, pac_n, "pac")
        ax_cy.set_title("Cytogenetic subtypes", fontsize=9.5, loc="left", pad=6,
                        fontweight="bold")
        ax_pa.set_title("PAC clusters", fontsize=9.5, loc="left", pad=6,
                        fontweight="bold")
        # combined size+colour legend (dot grows AND darkens with % positive)
        for v in (20, 50, 80):
            ax_leg.scatter([], [], s=s_of(v), color=cmap(norm(v)),
                           edgecolor="#3b5870", linewidth=0.4, label=f"{v}%")
        ax_leg.legend(title="% cells positive\n(≥1 raw read)", loc="center",
                      labelspacing=2.2, frameon=False, fontsize=8.5,
                      title_fontsize=8.5, borderpad=1.2, handletextpad=1.2)
        # (figure title/commentary in caption -- added by authors)
        save_fig(fig, HERE, "Figure7D_cyto_PAC_heatmap")
    print("\nDone -> Figure7D_cyto_PAC_heatmap.{pdf,png} + CSVs")


if __name__ == "__main__":
    main()
