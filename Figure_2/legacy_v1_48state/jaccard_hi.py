#!/usr/bin/env python
"""
High-resolution extension of the Jaccard stability sweep: push resolution up to
find where the number of stable clusters (per-cluster Jaccard >= 0.6) PLATEAUS.
Lean/RAM-safe: thread caps, explicit gc, B=20, high resolutions only.
Reuses the un-integrated atlas PCA (k=50, 50 PCs); igraph leiden n_iterations=2.
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
          "NUMBA_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[v] = "4"
import time, gc, warnings, numpy as np, pandas as pd, h5py, scanpy as sc, anndata as ad
from sklearn.metrics import adjusted_rand_score
warnings.filterwarnings("ignore"); sc.settings.verbosity = 0

OUT = "/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad/jaccard_hi"
os.makedirs(OUT, exist_ok=True)
_LOG = open(os.path.join(OUT, "run.log"), "a")
def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"; print(s, flush=True); _LOG.write(s+"\n"); _LOG.flush()

H5 = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ALSF_AML_plot/H5AD/ALSF_AML_Combo_3500_with_PAC_new.h5ad"
RES = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
B, FRAC, N_NEIGH, N_PCS = 20, 0.8, 50, 50
CHECKPOINT = 5
rng = np.random.default_rng(0)

f = h5py.File(H5, "r")
def col(name):
    node = f["obs"][name]
    if isinstance(node, h5py.Group):
        cats = [c.decode() if isinstance(c, bytes) else c for c in node["categories"][:]]
        return np.array([cats[i] if i >= 0 else "NA" for i in node["codes"][:]], dtype=object)
    arr = node[:]
    return np.array([a.decode() for a in arr], dtype=object) if arr.dtype.kind == "S" else arr
Xpca = f["obsm"]["X_pca"][:]
sub = col("AML Sub-Clusters"); sid = col("SampleID")
pac = col("Prognosis-Associated Clusters"); bc = col("_index")
f.close()

mask = sub != "NA"
Xl = np.ascontiguousarray(Xpca[mask][:, :N_PCS], dtype=np.float32)
subl, sidl, pacl, bcl = sub[mask], sid[mask], pac[mask], bc[mask]
n = Xl.shape[0]
log(f"leukemic cells {n}; HI Jaccard sweep B={B}, res={RES}, threads=4")

def leiden_multi(coords, res_list, seed):
    a = ad.AnnData(np.zeros((coords.shape[0], 1), dtype=np.float32)); a.obsm["X_pca"] = coords
    sc.pp.neighbors(a, n_neighbors=N_NEIGH, n_pcs=N_PCS, use_rep="X_pca", random_state=seed)
    out = {}
    for r in res_list:
        sc.tl.leiden(a, resolution=r, flavor="igraph", n_iterations=2, directed=False,
                     random_state=0, key_added="lx")
        _, out[r] = np.unique(a.obs["lx"].to_numpy(), return_inverse=True)
    del a; gc.collect()
    return out

t0 = time.time(); ref = leiden_multi(Xl, RES, seed=0)
K = {r: int(ref[r].max()+1) for r in RES}
log("reference %ds; nclust " % (time.time()-t0) + ", ".join(f"{r}:{K[r]}" for r in RES))
asg = pd.DataFrame({"barcode": bcl, "SampleID": sidl, "AML_SubCluster": subl, "PAC": pacl})
for r in RES: asg[f"res{r:.2f}"] = ref[r]
asg.to_csv(os.path.join(OUT, "reference_assignments.csv"), index=False)

jacc = {r: [[] for _ in range(K[r])] for r in RES}
ari  = {r: [] for r in RES}
pac_sub = {"AML_3":"FPAC_1","AML_24":"FPAC_2","AML_2":"PPAC_1","AML_11":"PPAC_2",
           "AML_14":"PPAC_3","AML_17":"PPAC_4","AML_23":"PPAC_5"}
for bi in range(B):
    idx = np.sort(rng.choice(n, size=int(FRAC*n), replace=False))
    ts = time.time(); lab = leiden_multi(Xl[idx], RES, seed=int(1000+bi))
    for r in RES:
        rl = ref[r][idx]; sl = lab[r]
        ari[r].append(adjusted_rand_score(rl, sl))
        kb = int(sl.max()+1)
        CM = np.zeros((K[r], kb), dtype=np.float64); np.add.at(CM, (rl, sl), 1.0)
        rs = CM.sum(1); cs = CM.sum(0)
        with np.errstate(divide="ignore", invalid="ignore"):
            J = CM / (rs[:, None] + cs[None, :] - CM)
        J[~np.isfinite(J)] = 0.0
        best = J.max(1)
        for c in range(K[r]): jacc[r][c].append(best[c])
    del lab; gc.collect()
    log(f"bootstrap {bi+1}/{B} in {time.time()-ts:.0f}s")
    if (bi+1) % CHECKPOINT == 0 or (bi+1) == B:
        prow, grow = [], []
        for r in RES:
            nst6 = nst7 = 0
            for c in range(K[r]):
                cm = ref[r] == c; mj = float(np.mean(jacc[r][c]))
                vcp = pd.Series(sidl[cm]).value_counts()
                pacs = pd.Series([pac_sub.get(x, "other") for x in subl[cm]])
                pv = pacs[pacs != "other"].value_counts()
                topPAC = f"{pv.index[0]}:{pv.iloc[0]/cm.sum():.0%}" if len(pv) else "-"
                prow.append([r, c, int(cm.sum()), mj, float(np.std(jacc[r][c])),
                             float(vcp.iloc[0]/cm.sum()), vcp.index[0], int(len(vcp)), topPAC])
                nst6 += mj >= 0.6; nst7 += mj >= 0.75
            grow.append([r, K[r], int(nst6), int(nst7),
                         float(np.mean([np.mean(jacc[r][c]) for c in range(K[r])])),
                         float(np.mean(ari[r]))])
        pd.DataFrame(prow, columns=["res","cluster","n_cells","mean_jaccard","sd_jaccard",
            "top_patient_frac","top_patient","n_patients","top_PAC"]
            ).to_csv(os.path.join(OUT, "jaccard_percluster.csv"), index=False)
        gg = pd.DataFrame(grow, columns=["res","n_clusters","n_stable_j60","n_stable_j75","mean_jaccard","mean_ARI"])
        gg.to_csv(os.path.join(OUT, "jaccard_global.csv"), index=False)
        log(f"[B={bi+1}] stable(J>=.6) " + ", ".join(f"{r}:{int(gg.loc[gg.res==r,'n_stable_j60'].iloc[0])}/{K[r]}" for r in RES))
log("DONE")
