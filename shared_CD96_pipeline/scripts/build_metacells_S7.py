#!/usr/bin/env python
# ====================================================================
# build_metacells_S7.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_metacells_S7.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Figure S7 panel A — CD96 vs comparator co-expression on AML METACELLS.

Per-cell gene-gene correlation in scRNA is confounded by dropout (technical zeros
track sequencing depth -> genes share zeros for non-biological reasons). We remove
this by aggregating leukemic cells into metacells: within each patient, k-means on
the 50-dim PCA partitions cells into ~75-cell transcriptionally coherent pools;
raw counts are summed per metacell and log-normalised (log1p of CP10k). CD96-vs-
target Pearson/Spearman are then computed across metacells.

Writes (into ../source_data):
  figS17_points_singlecell.csv   metacell points (marker, cd96, val, subtype)
  figS17_corr_r.csv              single_cell rows REPLACED with metacell r; bulk rows kept

Deterministic (k-means seed fixed); no fabricated data.
Run:  ../../.venv/bin/python build_metacells_S7.py
"""
import os
import h5py
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans

HERE = os.path.dirname(os.path.abspath(__file__))
RD = os.path.join(HERE, "..", "source_data")
H5AD = os.path.join(HERE, "..", "..", "ALSF_AML_plot", "H5AD",
                    "ALSF_AML_Combo_3500_with_PAC_new.h5ad")
# the 12 panel-E antigens (= the genes in Figure 7B/C)
GENES = ["CD96", "CD9", "SUCNR1", "IL1RAP", "FLT3", "CD33", "IL3RA", "CLEC12A", "CD7", "TNFRSF4", "ABCA7", "ITGAX"]
OTHERS = [g for g in GENES if g != "CD96"]
GLABEL = {"CD96": "CD96", "CD9": "CD9", "SUCNR1": "SUCNR1", "IL1RAP": "IL1RAP", "FLT3": "FLT3",
          "CD33": "CD33", "IL3RA": "CD123", "CLEC12A": "CLL-1", "CD7": "CD7", "TNFRSF4": "TNFRSF4",
          "ABCA7": "ABCA7", "ITGAX": "ITGAX"}
SC_LABEL = {"BCR/ABL": "BCR::ABL1", "RUNX1/RUNX1T1": "t(8;21)",
            "t(2;3)(p15;q26.2)": "t(2;3) MECOM", "t(7;14)(q21;q32)": "t(7;14)",
            "PML/RARA": "PML-RARA", "CN": "CN", "CBFB/MYH11": "inv(16)",
            "Tri(8)": "Trisomy 8", "Tri(15)": "Trisomy 15", "MLLr": "KMT2Ar",
            "del7q": "del(7q)", "MYB/GATA1": "MYB-GATA1",
            "NUP98/NSD1": "NUP98-NSD1", "Tri(8)/MLLr": "Tri8 + KMT2Ar"}
LEUK = [f"AML_{i}" for i in range(1, 28)]
TARGET = 75   # target cells per metacell


def _dec(a):
    return [x.decode() if isinstance(x, bytes) else x for x in a]


def _cat(f, k):
    return np.array(pd.Categorical.from_codes(
        f["obs"][k]["codes"][:].astype(int),
        categories=_dec(f["obs"][k]["categories"][:])).astype(object))


def main():
    print("loading h5ad ...", flush=True)
    import scipy.sparse as sp
    with h5py.File(H5AD, "r") as f:
        sub = _cat(f, "AML Sub-Clusters")
        cyto = _cat(f, "Cytogenetic")
        samp = _cat(f, "Sample")
        pca = f["obsm"]["X_pca"][:]
        rx = f["raw"]["X"]
        X = sp.csr_matrix((rx["data"][:], rx["indices"][:], rx["indptr"][:]),
                          shape=tuple(rx.attrs["shape"]))
        lib = np.asarray(X.sum(axis=1)).ravel().astype(np.float64)
        lib[lib == 0] = 1.0
        Xc = X.tocsc()
        rk = f["raw"]["var"].attrs.get("_index", b"_index")
        rk = rk.decode() if isinstance(rk, bytes) else rk
        sym = _dec(f["raw"]["var"][rk][:])
        rawg = {g: np.asarray(Xc[:, sym.index(g)].todense(), np.float64).ravel()
                for g in GENES}
    leuk = np.isin(sub, LEUK)
    mc = []
    for s in sorted(set(samp[leuk])):
        cells = np.where(leuk & (samp == s))[0]
        n = len(cells)
        if n == 0:
            continue
        k = max(1, int(round(n / TARGET)))
        labels = (np.zeros(n, dtype=int) if k == 1 else
                  KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(pca[cells]))
        st = SC_LABEL.get(cyto[cells][0], "Other")
        for lb in range(labels.max() + 1):
            cc = cells[labels == lb]
            if len(cc) == 0:
                continue
            slib = lib[cc].sum()
            row = {"subtype": st, "ncells": len(cc)}
            for g in GENES:
                row[g] = float(np.log1p(rawg[g][cc].sum() / slib * 1e4))
            mc.append(row)
    M = pd.DataFrame(mc)
    print(f"metacells: {len(M)} (median {int(M.ncells.median())} cells)", flush=True)

    rrows, pts = [], []
    for g in OTHERS:
        r = stats.pearsonr(M["CD96"], M[g])[0]
        rho = stats.spearmanr(M["CD96"], M[g])[0]
        rrows.append(dict(modality="single_cell", marker=g, label=GLABEL[g],
                          pearson=r, spearman=rho, n=len(M)))
        pts.append(pd.DataFrame({"marker": GLABEL[g], "cd96": M["CD96"].values,
                                 "val": M[g].values, "subtype": M["subtype"].values}))
    pd.concat(pts, ignore_index=True).to_csv(
        os.path.join(RD, "figS17_points_singlecell.csv"), index=False)
    cur = pd.read_csv(os.path.join(RD, "figS17_corr_r.csv"))
    pd.concat([pd.DataFrame(rrows), cur[cur.modality == "bulk"]],
              ignore_index=True).to_csv(os.path.join(RD, "figS17_corr_r.csv"), index=False)
    print(pd.DataFrame(rrows)[["label", "pearson", "spearman"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
