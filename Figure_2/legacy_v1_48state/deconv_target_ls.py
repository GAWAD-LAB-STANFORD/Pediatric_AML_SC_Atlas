#!/usr/bin/env python
"""
First-pass local deconvolution of the 30 leukemic states across TARGET bulk.
Signature = 30-state marker mean-CPM (built from leukemic cells); NNLS per patient
(no upload). Outputs per-patient fractions joined to survival + cytogenetic subtype,
plus TARGET cross-patient detection and cytogenetic-mean tables. Preliminary: a
normal-reference-augmented CIBERSORTx run is the confirmatory step.
"""
import os, numpy as np, pandas as pd
from scipy.optimize import nnls
SC = "/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
SIG = f"{SC}/signatures_all"; OUT = f"{SC}/target_ls"; os.makedirs(OUT, exist_ok=True)
PKG = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated"

# ---- 30-state signature (marker genes) with LS names ----
genes = np.load(f"{SIG}/genes.npy", allow_pickle=True).astype(str)
M = np.load(f"{SIG}/meancpm_res4.50.npy")                       # 63 x G
pc = pd.read_csv(f"{SC}/jaccard_hi/jaccard_percluster.csv"); pc = pc[pc.res == 4.5]
stab = pc[pc.mean_jaccard >= 0.6].copy().sort_values("n_cells", ascending=False)
raws = stab.cluster.astype(int).tolist()                        # stable cluster ids, LS_1..LS_30 order
ls_names = [f"LS_{i+1}" for i in range(len(raws))]
Ms = M[raws]                                                    # 30 x G
order = np.argsort(-Ms, 0); top = order[0]
tv = Ms[top, np.arange(Ms.shape[1])]; sv = Ms[order[1], np.arange(Ms.shape[1])]; fc = (tv+1)/(sv+1)
ch = []
for c in range(len(raws)):
    gi = np.where((top == c) & (tv > 1))[0]
    if len(gi): ch.append(gi[np.argsort(-fc[gi])][:50])
sig_idx = np.unique(np.concatenate(ch))
S_genes = genes[sig_idx]; S = Ms[:, sig_idx].T                  # (sig_genes x 30) linear CPM

# ---- TARGET mixture ----
mix = pd.read_csv(f"{PKG}/cibersortx_input/TARGET_mixture_symbols_CPM.txt", sep="\t", index_col=0)
common = [g for g in S_genes if g in mix.index]
print(f"signature genes {len(S_genes)}, in TARGET {len(common)}", flush=True)
Sc = pd.DataFrame(S, index=S_genes, columns=ls_names).loc[common].values   # (n x 30)
Y = mix.loc[common].values                                     # (n x patients)
pts = mix.columns.tolist()

# ---- NNLS per patient ----
F = np.zeros((len(pts), len(raws)))
for j in range(len(pts)):
    x, _ = nnls(Sc, Y[:, j]); s = x.sum(); F[j] = x/s if s > 0 else x
frac = pd.DataFrame(F, columns=ls_names); frac.insert(0, "sample", pts)
frac.to_csv(f"{OUT}/LS_target_fractions.csv", index=False)

# ---- join survival + cytogenetic subtype ----
surv = pd.read_csv(f"{PKG}/__SUBMISSION_PACKAGE/08_Source_Data/Figure_2/Figure2E_TARGET_survival_input.csv")
cd = pd.read_csv(f"{PKG}/__SUBMISSION_PACKAGE/08_Source_Data/Figure_2/Figure2CD_TARGET_PAC_perpatient.csv")[["sample","subtype2"]]
J = frac.merge(surv[["sample","os_time","os_event","efs_time","efs_event","age","wbc","risk"]], on="sample", how="left")
J = J.merge(cd.drop_duplicates("sample"), on="sample", how="left")
J.to_csv(f"{OUT}/LS_target_survival_input.csv", index=False)
print(f"patients: {len(J)}, with survival: {J.os_time.notna().sum()}", flush=True)

# ---- TARGET cross-patient detection ----
det = []
for ls in ls_names:
    v = frac[ls].values
    det.append([ls, round((v > 0.01).mean(), 3), round((v > 0.05).mean(), 3), round(float(v.mean()), 4), round(float(np.median(v)), 4)])
det = pd.DataFrame(det, columns=["LS","frac_pts>1%","frac_pts>5%","mean_frac","median_frac"])
det.to_csv(f"{OUT}/LS_target_detection.csv", index=False)

# ---- cytogenetic mean fraction + sanity (argmax subtype per state) ----
jc = J.dropna(subset=["subtype2"])
cyto_mean = jc.groupby("subtype2")[ls_names].mean().T
cyto_mean.to_csv(f"{OUT}/LS_target_cyto_meanfrac.csv")
print("\n=== TARGET detection (top/bottom by cross-patient presence) ===")
print(det.sort_values("frac_pts>1%", ascending=False).to_string(index=False))
print("\n=== sanity: TARGET subtype where each state peaks (expect match to discovery cytogenetics) ===")
for ls in ls_names:
    row = cyto_mean.loc[ls]; print(f"  {ls}: peaks in {row.idxmax()} ({row.max():.3f})")
print("\nDONE -> ", OUT)
