#!/usr/bin/env python
# ====================================================================
# build_FigureS_CD96_correlations.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_FigureS_CD96_correlations.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Supplement: does CD96 co-express with the other 8 surface targets?

For each of the 8 comparators, CD96 (x) vs that target (y), in TWO modalities:
  rows 1-2  -- single AML cells (scRNA leukemic cells AML_1..27; CP10K+log1p)
  rows 3-4  -- TARGET bulk primary-AML samples (log2 CPM)
Points coloured by cytogenetic subtype; a Pearson r (CD96 vs target) is printed
on every panel (Pearson + Spearman + n saved to CSV).

Interpretation: low / negative CD96-vs-comparator correlation = CD96 marks
DIFFERENT cells/patients than the clinical targets (complementary), not
redundant.  Caveats: single-cell r is inflated by shared dropout (co-zeros);
bulk uses all primary AML so normal-cell content can add shared variance.
Pure data (h5ad + TARGET STAR counts); no fabricated values.
"""
import os
import sys
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, "supplements")
sys.path.insert(0, SUPP)
from _fig7_style import PALETTE, RC, save_fig

H5AD = os.path.join(HERE, "..", "ALSF_AML_plot", "H5AD",
                    "ALSF_AML_Combo_3500_with_PAC_new.h5ad")
FILT = os.path.join(HERE, "..", "ALSF_AML_plot", "H5AD",
                    "NCBI_TARGET_Star_counts_filtered.csv")
MERGED = os.path.join(HERE, "Figure7F_TARGET_CD96_merged.csv")
BLAST = "Bone marrow leukemic blast percentage (%)"

ENS = {"CD96": "ENSG00000153283", "CD33": "ENSG00000105383",
       "IL3RA": "ENSG00000185291", "CLEC12A": "ENSG00000172322",
       "FLT3": "ENSG00000122025", "MSLN": "ENSG00000102854",
       "CD70": "ENSG00000125726", "ADGRG1": "ENSG00000205336",
       "IL1RAP": "ENSG00000196083"}
GENES = list(ENS)
OTHERS = [g for g in GENES if g != "CD96"]          # 8 comparators
GLABEL = {"CD96": "CD96", "CD33": "CD33", "IL3RA": "CD123 (IL3RA)",
          "CLEC12A": "CLL-1 (CLEC12A)", "FLT3": "FLT3", "MSLN": "MSLN",
          "CD70": "CD70", "ADGRG1": "GPR56 (ADGRG1)", "IL1RAP": "IL1RAP"}
LEUK_SUBS = [f"AML_{i}" for i in range(1, 28)]

SUBORDER = ["CN", "t(8;21)", "inv(16)", "KMT2Ar", "PML-RARA", "NUP98r", "Other"]
SUBCOL = {"CN": "#0072B2", "t(8;21)": "#009E73", "inv(16)": "#E69F00",
          "KMT2Ar": "#c0392b", "PML-RARA": "#CC79A7", "NUP98r": "#56B4E9",
          "Other": "#b9c2cb"}
SC_MAP = {"CN": "CN", "RUNX1/RUNX1T1": "t(8;21)", "CBFB/MYH11": "inv(16)",
          "MLLr": "KMT2Ar", "PML/RARA": "PML-RARA", "NUP98/NSD1": "NUP98r"}


def bulk_sub(row):
    fz = str(row.get("Gene Fusion", "")).strip()
    if fz and fz.lower() != "nan":
        if fz.startswith("KMT2A-") or fz == "other MLL":
            return "KMT2Ar"
        if fz == "RUNX1-RUNX1T1":
            return "t(8;21)"
        if fz == "CBFB-MYH11":
            return "inv(16)"
        if fz.startswith("NUP98-"):
            return "NUP98r"
        return "Other"
    return {"Normal": "CN", "t(8;21)": "t(8;21)", "inv(16)": "inv(16)",
            "PML-RARA": "PML-RARA", "MLL": "KMT2Ar"}.get(
        str(row.get("Primary Cytogenetic Code", "")).strip(), "Other")


def _decode(a):
    return [x.decode() if isinstance(x, bytes) else x for x in a]


def _cat(f, k):
    return np.array(pd.Categorical.from_codes(
        f["obs"][k]["codes"][:].astype(int),
        categories=_decode(f["obs"][k]["categories"][:])).astype(object))


def load_singlecell():
    with h5py.File(H5AD, "r") as f:
        cyto = _cat(f, "Cytogenetic")
        sub = _cat(f, "AML Sub-Clusters")
        rx = f["raw"]["X"]
        X = sp.csr_matrix((rx["data"][:], rx["indices"][:], rx["indptr"][:]),
                          shape=tuple(rx.attrs["shape"]))
        rk = f["raw"]["var"].attrs.get("_index", b"_index")
        rk = rk.decode() if isinstance(rk, bytes) else rk
        symbols = _decode(f["raw"]["var"][rk][:])
        lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float32)
        lib[lib == 0] = 1.0
        Xc = X.tocsc()
        cols = {}
        for g in GENES:
            v = np.asarray(Xc[:, symbols.index(g)].todense(), np.float32).ravel()
            cols[g] = np.log1p(v / lib * 1e4)              # CP10K + log1p
    d = pd.DataFrame(cols)
    d["sub"] = sub
    d["subtype"] = pd.Series(cyto).map(SC_MAP).fillna("Other").values
    return d[d["sub"].isin(LEUK_SUBS)].reset_index(drop=True)


def load_bulk():
    raw = pd.read_csv(FILT).set_index("ensembl_gene_id")
    raw.index = raw.index.astype(str).str.split(".").str[0]
    g9 = raw.loc[[ENS[g] for g in GENES]].T
    g9.columns = GENES
    m = pd.read_csv(MERGED).set_index("sample")
    ids = [s for s in g9.index if s in m.index and pd.notna(m.loc[s, "lib_size"])]
    lib = m.loc[ids, "lib_size"].astype(float)
    d = np.log2(g9.loc[ids].div(lib, axis=0) * 1e6 + 1).reset_index(drop=True)
    d["subtype"] = [bulk_sub(m.loc[s]) for s in ids]
    return d


def scatter(ax, df, marker, unit, max_pts=20000):
    x = df["CD96"].to_numpy(float); y = df[marker].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, sub = x[ok], y[ok], df["subtype"].to_numpy()[ok]
    pear = stats.pearsonr(x, y)[0] if len(x) > 2 else np.nan
    spear = stats.spearmanr(x, y)[0] if len(x) > 2 else np.nan
    # plot (downsample for rendering only; r computed on all points)
    idx = np.arange(len(x))
    if len(idx) > max_pts:
        idx = np.linspace(0, len(x) - 1, max_pts).astype(int)
    for s in SUBORDER:                                     # Other first, hot last
        m = sub[idx] == s
        if m.any():
            ax.scatter(x[idx][m], y[idx][m], s=5, c=SUBCOL[s],
                       alpha=0.35, linewidths=0, zorder=2 if s != "Other" else 1)
    rc = PALETTE["CD96"] if abs(pear) < 0.3 else "#1c2833"
    ax.text(0.04, 0.96, f"r = {pear:+.2f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=9, fontweight="bold", color=rc,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))
    ax.set_xlabel(f"CD96 {unit}", fontsize=7.5)
    ax.set_ylabel(f"{GLABEL[marker]} {unit}", fontsize=7.5)
    ax.tick_params(labelsize=6.5)
    for sp_ in ("top", "right"):
        ax.spines[sp_].set_visible(False)
    return dict(marker=marker, pearson=pear, spearman=spear, n=int(len(x)))


def main():
    print("  loading single cells ..."); sc = load_singlecell()
    print(f"   leukemic cells: {len(sc):,}")
    print("  loading bulk ...");          bk = load_bulk()
    print(f"   TARGET samples: {len(bk):,}")

    rows = []
    with mpl.rc_context(RC):
        fig, axes = plt.subplots(4, 4, figsize=(15.0, 15.2))
        for i, g in enumerate(OTHERS):                    # single cell rows 0-1
            r = scatter(axes[i // 4, i % 4], sc, g, "(log-norm)")
            r["modality"] = "single_cell"; rows.append(r)
        for i, g in enumerate(OTHERS):                    # bulk rows 2-3
            r = scatter(axes[2 + i // 4, i % 4], bk, g, "(log2 CPM)")
            r["modality"] = "bulk"; rows.append(r)

        # block labels
        fig.text(0.012, 0.985, "A   Single AML cells  (scRNA, leukemic cells "
                 f"n={len(sc):,})", fontsize=12, fontweight="bold",
                 color="#1c2833", ha="left", va="top")
        fig.text(0.012, 0.495, "B   TARGET bulk  (primary AML "
                 f"n={len(bk):,})", fontsize=12, fontweight="bold",
                 color="#1c2833", ha="left", va="top")
        # subtype legend
        handles = [plt.Line2D([], [], marker="o", ls="", ms=7,
                              mfc=SUBCOL[s], mec="none", label=s)
                   for s in SUBORDER]
        fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False,
                   fontsize=9, bbox_to_anchor=(0.5, -0.005),
                   title="cytogenetic subtype", title_fontsize=9)
        # (figure title/commentary in caption -- added by authors)
        fig.subplots_adjust(left=0.06, right=0.99, top=0.955, bottom=0.055,
                            hspace=0.42, wspace=0.32)
        save_fig(fig, HERE, "FigureS_CD96_correlations")

    out = pd.DataFrame(rows)[["modality", "marker", "pearson", "spearman", "n"]]
    out.round(3).to_csv(os.path.join(HERE, "FigureS_CD96_correlations.csv"),
                        index=False)
    print("\n", out.round(2).to_string(index=False))
    print("\nDone -> FigureS_CD96_correlations.{pdf,png,csv}")


if __name__ == "__main__":
    main()
