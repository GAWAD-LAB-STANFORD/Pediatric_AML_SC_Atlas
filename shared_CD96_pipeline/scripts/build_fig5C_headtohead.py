#!/usr/bin/env python
"""Figure 5 panel C (head-to-head, corrected scheme): efficacy vs toxicity for the surfaceome candidates
overlaid with the clinically-pursued AML antigens + FLT3 control. Recomputed from the MAIN atlas raw
counts (h5py) on the corrected leukemic definition (Cell_Type in config.LEUK_CELLTYPES AND SampleType==
"AML") with normal HSPC/Myeloid_Pro taken from Cell_Type, so panel C is consistent with panels A/D.
For every gene: x = leak (max % normal HSPC/Myeloid_Pro positive), y = % leukemic cells positive,
size = % of the 28 AML patients with >=20% leukemic positivity. Composite (candidate ranking) = mean of
% leukemic + % patients (2-axis; the retired PPAC axis removed). No fabricated values."""
import os, sys, numpy as np, pandas as pd, scipy.sparse as sp, h5py
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C

SURF = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ref/surfaceome_genes.txt"
OUT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/shared_CD96_pipeline/source_data/fig5C_headtohead.csv"
GATE = 10.0
LEADS = ["CD96", "CD9", "SUCNR1", "IL1RAP"]
CLIN6 = {"CD33": "CD33", "IL3RA": "CD123", "CLEC12A": "CLL-1", "CD38": "CD38", "CD7": "CD7", "CD70": "CD70"}
CONTROL = {"FLT3": "FLT3"}
# former "PPAC-specific" markers -> now the "additional candidate" overlay group (label kept as 'ppac' for
# fig5.R which recolours ppac->additional/candidate and CD7->clinical).
PPAC = {"CD7": "CD7", "TNFRSF4": "TNFRSF4", "ABCA7": "ABCA7", "ITGAX": "ITGAX"}
TCR = {"CD3D", "CD3E", "CD3G", "CD247", "TRAC", "TRBC1", "TRBC2", "TRGC1", "TRGC2", "TRDC"}
TOPN_CANDIDATES = 40

dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
def cat(f, k):
    o = f["obs"][k]; return dec(o["categories"][:])[o["codes"][:].astype(int)]

f = h5py.File(C.H5, "r")
ct = cat(f, "Cell Type"); st = cat(f, "SampleType"); sid = cat(f, "SampleID")
rk = f["raw"]["var"].attrs.get("_index", b"_index"); rk = rk.decode() if isinstance(rk, bytes) else rk
sym = list(dec(f["raw"]["var"][rk][:]))
rx = f["raw"]["X"]
X = sp.csr_matrix((rx["data"][:], rx["indices"][:], rx["indptr"][:]), shape=tuple(rx.attrs["shape"]))
f.close()

surf = [g.strip() for g in open(SURF) if g.strip()]
need = sorted(set([g for g in surf if g in sym]) | set(CLIN6) | set(CONTROL) | set(PPAC))
present = [g for g in need if g in sym]
gi = [sym.index(g) for g in present]
posS = (X[:, gi] > 0).tocsc()                                  # cells x genes, boolean
idx = {g: k for k, g in enumerate(present)}

leuk = np.isin(ct, C.LEUK_CELLTYPES) & (st == "AML")           # corrected leukemic set
hspc = ct == "0_HSPC"; mye = ct == "Myeloid_Pro"              # normal marrow reference
patients = sorted(set(sid[leuk])); N_PT = len(patients)
hspc_pct = 100 * np.asarray(posS[hspc].mean(0)).ravel()
mye_pct  = 100 * np.asarray(posS[mye].mean(0)).ravel()
aml_pct  = 100 * np.asarray(posS[leuk].mean(0)).ravel()
npts = np.zeros(len(present), int)
for s in patients:
    npts += (np.asarray(posS[leuk & (sid == s)].mean(0)).ravel() >= 0.20)

def metrics(i):
    comp = (aml_pct[i] + npts[i] / N_PT * 100) / 2             # 2-axis composite (no PPAC term)
    return dict(leak=round(max(hspc_pct[i], mye_pct[i]), 2), scRNA_AML_pct=round(aml_pct[i], 2),
                scRNA_HSPC_pct=round(hspc_pct[i], 2), scRNA_Myeloid_pct=round(mye_pct[i], 2),
                n_pts_pct=round(npts[i] / N_PT * 100, 1), composite=round(comp, 2))

rows = {}
gated = []
for g in present:
    if g in TCR: continue
    mm = metrics(idx[g])
    if mm["scRNA_HSPC_pct"] <= GATE and mm["scRNA_Myeloid_pct"] <= GATE:
        gated.append((g, mm))
gated.sort(key=lambda t: -t[1]["composite"])
for rank, (g, mm) in enumerate(gated[:TOPN_CANDIDATES], 1):
    rows[g] = dict(gene=g, label=g, **mm, klass=("lead" if g in LEADS else "candidate"), comp_rank=rank)
for g, lab in CONTROL.items():
    rows[g] = dict(gene=g, label=lab, **metrics(idx[g]), klass="control", comp_rank=999)
for g, lab in CLIN6.items():
    prev = rows.get(g); cr = prev["comp_rank"] if prev else 999
    rows[g] = dict(gene=g, label=lab, **metrics(idx[g]),
                   klass=("lead" if g in LEADS else "clinical"), comp_rank=cr)
for g, lab in PPAC.items():                                    # overlay last (incl. CD7)
    prev = rows.get(g); cr = prev["comp_rank"] if prev else 999
    rows[g] = dict(gene=g, label=lab, **metrics(idx[g]), klass="ppac", comp_rank=cr)

df = pd.DataFrame(rows.values())
df["do_label"] = (df.klass.isin(["lead", "clinical", "control", "ppac"])) | (df.comp_rank <= 10)
df = df.sort_values(["klass", "composite"], ascending=[True, False])
for d in [OUT, OUT.replace("07_Code/shared_CD96_pipeline/source_data", "08_Source_Data/Figure_5")]:
    df.to_csv(d, index=False)
print("wrote fig5C_headtohead.csv (corrected: %d AML patients, leukemic=%s & AML)" % (N_PT, C.LEUK_CELLTYPES))
print(df[df.klass.isin(['clinical', 'control', 'lead', 'ppac'])][
    ['label', 'klass', 'leak', 'scRNA_AML_pct', 'n_pts_pct']].to_string(index=False))
print("\nCD96:", metrics(idx["CD96"]))
