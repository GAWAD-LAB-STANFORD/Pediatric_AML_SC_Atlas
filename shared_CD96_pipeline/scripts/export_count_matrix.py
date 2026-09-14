#!/usr/bin/env python
# ====================================================================
# export_count_matrix.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python export_count_matrix.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Export the RAW COUNT matrix from the master h5ad into the package as a compact,
counts-only AnnData (+ cell-metadata CSV). This is the primary raw input to
reproduce the scRNA figures; the 44 GB master object's derived layers (scaled X,
imputed, log-norm) are NOT included — they are regenerable from these counts.

Source : layers['Raw_Counts']  (integer counts, 96,627 cells x 27,346 genes)
Output : count_matrix/ALSF_AML_raw_counts.h5ad  + count_matrix/cell_metadata.csv.gz
Includes obs (all cell annotations) + obsm X_umap/X_pca so the object is directly
usable (UMAP plots, neighbors, metacells) without recomputation.
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import anndata as ad

HERE = os.path.dirname(os.path.abspath(__file__))
H5AD = os.path.join(HERE, "..", "..", "ALSF_AML_plot", "H5AD",
                    "ALSF_AML_Combo_3500_with_PAC_new.h5ad")
OUT = os.path.join(HERE, "..", "count_matrix")
os.makedirs(OUT, exist_ok=True)


def dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else x for x in a])


def read_col(o):
    if isinstance(o, h5py.Group) and "categories" in o:           # categorical
        return pd.Categorical.from_codes(o["codes"][:].astype(int), dec(o["categories"][:]))
    arr = o[:]
    return dec(arr) if (arr.dtype.kind in ("S", "O")) else arr


with h5py.File(H5AD, "r") as f:
    rc = f["layers"]["Raw_Counts"]
    X = sp.csr_matrix((rc["data"][:], rc["indices"][:], rc["indptr"][:]),
                      shape=tuple(rc.attrs["shape"]))
    og = f["obs"]
    ik = og.attrs.get("_index", b"_index"); ik = ik.decode() if isinstance(ik, bytes) else ik
    cells = dec(og[ik][:])
    obs = {}
    for k in og.keys():
        if k == ik:
            continue
        try:
            col = read_col(og[k])
            if len(col) == len(cells):
                obs[k] = col
        except Exception:
            pass
    vk = f["var"].attrs.get("_index", b"_index"); vk = vk.decode() if isinstance(vk, bytes) else vk
    genes = dec(f["var"][vk][:])
    obsm = {k: f["obsm"][k][:] for k in ("X_umap", "X_pca") if "obsm" in f and k in f["obsm"]}

obs_df = pd.DataFrame(obs, index=pd.Index(cells, name="cell_id"))
A = ad.AnnData(X=X, obs=obs_df, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
for k, v in obsm.items():
    A.obsm[k] = v
print("AnnData:", A.shape, "| nnz", X.nnz, "| obs cols", obs_df.shape[1], "| obsm", list(obsm))

A.write_h5ad(os.path.join(OUT, "ALSF_AML_raw_counts.h5ad"), compression="gzip")
obs_df.to_csv(os.path.join(OUT, "cell_metadata.csv.gz"))
print("wrote ->", OUT)
