#!/usr/bin/env python
"""Figure 2B (v3, corrected data): leukemic states coloured by whole-cohort EFS
prognosis on the full-atlas UMAP. Leukemic cells = Cell Type in
{AML,AML-MKI67,AML-PCNA,AML-CD1C} (the corrected 70,108-cell set), assigned to the
res-4.5 stable states (Jaccard>=0.35) from reference_assignments.csv; each state is
coloured favorable / poor (whole-cohort EFS Cox FDR<0.10) or grey (ns / non-stable).
Prognosis is EFS-based (unified). No fabricated values."""
import os, sys, numpy as np, pandas as pd, h5py, anndata as ad, scanpy as sc, matplotlib.pyplot as plt
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C
D2 = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
OUT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
os.makedirs(OUT, exist_ok=True)

# --- state -> prognosis (whole-cohort EFS) ---
prog = pd.read_csv(os.path.join(D2, "LS_wholecohort_prognosis.csv"))
ls2prog = dict(zip(prog.LS, prog.prognosis))          # favorable / poor / ns

# --- per-cell stable-state labels from the corrected reference clustering ---
ra = pd.read_csv(C.REFA)
ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, lslabs = C.stable_raws()                         # res4.50 cluster -> LS_i
raw2ls = {rw: lslabs[i] for i, rw in enumerate(raws)}
ra["LS"] = ra["res4.50"].map(raw2ls)                   # None for non-stable clusters
ra["prog"] = ra["LS"].map(ls2prog).fillna("ns")        # ns for non-stable / ns states

# sanity: healthy-BM contamination + stable coverage
hbm = ra["SampleID"].astype(str).str.startswith("0_Healthy").mean()
cov = ra["LS"].notna().mean()
print(f"leukemic cells: {len(ra)}; healthy-BM fraction: {hbm:.3f}; in a stable state: {cov:.3f}")
print(f"cells favorable {int((ra.prog=='favorable').sum())} / poor {int((ra.prog=='poor').sum())} / ns {int((ra.prog=='ns').sum())}")

# --- UMAP coords + Cell Type for ALL atlas cells; match leukemic by barcode ---
f = h5py.File(C.H5, "r")
obs_names = [x.decode() if isinstance(x, bytes) else x for x in f["obs"]["_index"][:]]
umap = f["obsm/X_umap"][:]
ctn = f["obs"]["Cell Type"]                                   # normal-HSC reference = Cell Type "0_HSPC"
ctcats = [c.decode() if isinstance(c, bytes) else c for c in ctn["categories"][:]]
celltype = np.array([ctcats[i] if i >= 0 else "NA" for i in ctn["codes"][:]], dtype=object)
f.close()
N = len(obs_names)
pos = {b: i for i, b in enumerate(obs_names)}
cellprog = np.array(["background"] * N, dtype=object)
miss = 0
for bc, pr in zip(ra["barcode"], ra["prog"]):
    i = pos.get(bc)
    if i is None: miss += 1; continue
    cellprog[i] = pr
print(f"barcodes not matched to h5ad: {miss}")
# normal-HSC reference overlay: MARKER-GATED bona-fide HSC (HLF/AVP+, CD34+, negative for
# LYZ/CD14/GATA1/MS4A1) -- the atlas "0_HSPC" label is a broad progenitor cluster with ~2/3
# lineage-marker+ cells, so we gate to the primitive HSC fraction. LANDMARK, not a prognostic
# category (NORM_HSPC is not prognostic). Barcodes from normal_HSC_gated_barcodes.csv.
gated = set(pd.read_csv(os.path.join(D2, "normal_HSC_gated_barcodes.csv"))["barcode"].astype(str))
cellprog[np.array([b in gated for b in obs_names])] = "normal HSC"
print(f"normal HSC (marker-gated) cells overlaid: {int((cellprog=='normal HSC').sum())}")

# --- plot: favorable teal, poor red, other leukemic (ns) grey, normal HSC orange, background faint ---
order = ["background", "ns", "favorable", "poor", "normal HSC"]  # draw normal-HSC landmark on top
colmap = {"background": "#eeeeee", "ns": "#c9c9c9", "favorable": "#0A9396", "poor": "#AE2012", "normal HSC": "#EE9B00"}
fig, axx = plt.subplots(figsize=(6.2, 5.6))
for k in order:
    m = cellprog == k
    axx.scatter(umap[m, 0], umap[m, 1], s=(1.2 if k in ("background","ns") else (9.0 if k=="normal HSC" else 3.0)),
                c=colmap[k], linewidths=0, rasterized=True,
                label={"background":"non-leukemic","ns":"leukemic (n.s.)",
                       "favorable":"favorable EFS","poor":"poor EFS","normal HSC":"normal HSC (reference)"}[k],
                alpha=(0.35 if k=="background" else 0.95))
axx.set_xticks([]); axx.set_yticks([]); axx.set_frame_on(False)
lg = axx.legend(loc="lower left", frameon=False, markerscale=3, fontsize=9, handletextpad=0.2)
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"Figure_2B_prognosis_umap.{ext}"), dpi=200, bbox_inches="tight")
plt.close()
print("wrote Figure_2B_prognosis_umap.{png,pdf} ->", OUT)
