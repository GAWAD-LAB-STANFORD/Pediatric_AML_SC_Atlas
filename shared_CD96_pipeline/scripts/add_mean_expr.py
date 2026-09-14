#!/usr/bin/env python
# ====================================================================
# add_mean_expr.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python add_mean_expr.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Recompute the dot-plot FILL channel (mean expression) for the scRNA dot plots
and merge it into the source-data CSVs. Closes the reproducibility gap flagged in
the audit: the upstream builders emitted % positive only, and the dot FILL = mean
expression had no committed builder.

mean_expr = mean over all cells in the group of log1p(raw / total_counts * 1e4)
            (standard scanpy dot-plot colour; zeros included).
pct       = 100 * fraction of cells in the group with raw count > 0.
No randomness, no fabricated values — pure recompute from h5ad raw counts.

Recipes recovered by matching the deployed `pct` column exactly:
  * fig7B_spec.csv  (Fig 7A specificity): group = compartment (Cell Type ->
      CELLTYPE_GROUPING), all cells.  REPRODUCES THE DEPLOYED pct EXACTLY
      (max|diff| < 0.01) -> this file is rewritten in place.
  * fig5D_pac.csv  (Fig 6D PAC):  group = obs 'PAC_anno', all cells.
  * fig5D_cyto.csv (Fig 6D cyto): group = obs 'Cytogenetic', core-AML cells
      (Cell Type in {AML, AML-MKI67, AML-PCNA}).
  fig6D_pac/cyto reproduce the deployed pct to ~1-3% for the multi-patient
  karyotypes (a per-sample/downsampling nuance of the original, now-drifted
  builder) and exactly for the rest; their mean_expr is therefore reported but
  NOT overwritten by default (pass --write-6d to regenerate them from this clean
  recipe, which shifts a few 6D pct values by <=3%).
"""
import os, sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
H5AD = os.path.join(HERE, "..", "..", "ALSF_AML_plot", "H5AD",
                    "ALSF_AML_Combo_3500_with_PAC_new.h5ad")
SRC = os.path.join(HERE, "..", "source_data")
GENES = ["CD96", "CD33", "IL3RA", "CLEC12A", "FLT3", "MSLN", "CD70", "ADGRG1", "IL1RAP"]
AMLCORE = {"AML", "AML-MKI67", "AML-PCNA"}
CELLTYPE_GROUPING = {
    "AML": "AML (leukemic)", "AML-MKI67": "AML (leukemic)", "AML-PCNA": "AML (leukemic)",
    "0_HSPC": "HSPC", "Myeloid_Pro": "Myeloid_Pro",
    "CTL": "T cells", "Naïve_CD4T": "T cells", "Naïve_CD8T": "T cells",
    "Activated_CD4T": "T cells", "AML-CTL": "T cells", "AML-CD4T": "T cells",
    "AML-Naïve_CD8T": "T cells", "NK": "NK", "AML-NK": "NK",
    "CD34+ProB": "B-lineage", "ProB": "B-lineage", "PreB": "B-lineage",
    "CD20+B": "B-lineage", "PlasmaB": "B-lineage", "AML-B": "B-lineage",
    "CD14_Monocyte": "Mature myeloid", "CD16_Monocyte": "Mature myeloid",
    "Macrophage": "Mature myeloid", "mDC": "Mature myeloid", "pDC": "Mature myeloid",
    "AML-CD14": "Mature myeloid", "AML-CD1C": "Mature myeloid",
    "Erythrocytes": "Erythroid", "AML-Ery": "Erythroid",
}


def _decode(a):
    return [x.decode() if isinstance(x, bytes) else x for x in a]


def _cat(f, k):
    o = f["obs"][k]
    c = np.array(_decode(o["categories"][:]))
    return c[o["codes"][:].astype(int)]


def main(write_6d=False):
    print("loading h5ad ...")
    with h5py.File(H5AD, "r") as f:
        ctype = _cat(f, "Cell Type"); cyto = _cat(f, "Cytogenetic")
        pacA = _cat(f, "PAC_anno")
        rx = f["raw"]["X"]
        X = sp.csr_matrix((rx["data"][:], rx["indices"][:], rx["indptr"][:]),
                          shape=tuple(rx.attrs["shape"]))
        rk = f["raw"]["var"].attrs.get("_index", b"_index")
        rk = rk.decode() if isinstance(rk, bytes) else rk
        symbols = list(_decode(f["raw"]["var"][rk][:]))
    total = np.asarray(X.sum(axis=1)).ravel().astype(np.float64); total[total == 0] = 1.0
    comp = np.array([CELLTYPE_GROUPING.get(c) for c in ctype], dtype=object)
    amlcore = np.isin(ctype, list(AMLCORE))
    Xc = X.tocsc(); expr, pos = {}, {}
    for g in GENES:
        raw = np.asarray(Xc[:, symbols.index(g)].todense(), np.float64).ravel()
        expr[g] = np.log1p(raw / total * 1e4); pos[g] = raw > 0

    def process(fname, gcol, cells, mask, tol, write):
        path = os.path.join(SRC, fname); df = pd.read_csv(path)
        me, pc = [], []
        for _, r in df.iterrows():
            m = (cells == r[gcol]) & mask
            me.append(round(float(expr[r["gene"]][m].mean()), 5) if m.sum() else np.nan)
            pc.append(float(100.0 * pos[r["gene"]][m].mean()) if m.sum() else np.nan)
        pdiff = np.nanmax(np.abs(np.array(pc) - df["pct"].values))
        mdiff = np.nanmax(np.abs(np.array(me) - df["mean_expr"].values)) if "mean_expr" in df else np.nan
        print(f"  {fname}: max|pct recompute - deployed|={pdiff:.3f}  max|mean_expr diff|={mdiff:.4f}")
        if write:
            assert pdiff < tol, f"{fname} pct mismatch {pdiff} > tol {tol}"
            df["mean_expr"] = me; df.to_csv(path, index=False)
            print(f"    -> rewrote {fname} (mean_expr reproduced from h5ad)")
        else:
            print(f"    -> left {fname} unchanged (deployed values retained; --write-6d to regenerate)")

    allmask = np.ones(X.shape[0], bool)
    process("fig7B_spec.csv", "compartment", comp, allmask, 0.01, write=True)
    process("fig5D_pac.csv", "group", pacA, allmask, 5.0, write=write_6d)
    process("fig5D_cyto.csv", "group", cyto, amlcore, 5.0, write=write_6d)
    print("done.")


if __name__ == "__main__":
    main(write_6d=("--write-6d" in sys.argv))
