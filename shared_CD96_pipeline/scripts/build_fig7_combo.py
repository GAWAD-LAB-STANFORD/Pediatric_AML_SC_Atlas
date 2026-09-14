#!/usr/bin/env python
"""Figure 7 — additive efficacy vs marrow toxicity of CD96 + one partner. Corrected scheme: recomputed
from the MAIN atlas raw counts (h5py) on the corrected leukemic definition (Cell_Type in
config.LEUK_CELLTYPES AND SampleType=="AML") with normal marrow = Cell_Type 0_HSPC/Myeloid_Pro.
WITHIN patients : % leukemic cells covered by (CD96 OR partner), mean across CD96-targetable patients.
BETWEEN patients: % patients with >=20% of blasts covered by (CD96 OR partner).
TOXICITY (x)    : % of normal HSPC/Myeloid_Pro cells positive for (CD96 OR partner).
-> source_data/fig7_combo.csv (+ 08_Source_Data mirror). No fabricated values."""
import os, sys, numpy as np, pandas as pd, scipy.sparse as sp, h5py
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C

SD = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/shared_CD96_pipeline/source_data"
SD08 = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_5"
GENES = {"CD96": "CD96", "CD9": "CD9", "SUCNR1": "SUCNR1", "IL1RAP": "IL1RAP", "FLT3": "FLT3",
         "CD33": "CD33", "IL3RA": "CD123", "CLEC12A": "CLL-1", "CD7": "CD7", "TNFRSF4": "TNFRSF4",
         "ABCA7": "ABCA7", "ITGAX": "ITGAX"}
GROUP = {"CD96": "lead", "CD9": "lead", "SUCNR1": "lead", "IL1RAP": "lead", "FLT3": "control",
         "CD33": "clinical", "CD123": "clinical", "CLL-1": "clinical", "CD7": "ppac",
         "TNFRSF4": "ppac", "ABCA7": "ppac", "ITGAX": "ppac"}
HIGHLY_TOXIC = {"CD9", "CD123", "IL1RAP", "ABCA7"}    # organ-toxic in Fig 6 -> red ring
COVER, CD96POS = 20.0, 10.0

dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
def cat(f, k):
    o = f["obs"][k]; return dec(o["categories"][:])[o["codes"][:].astype(int)]

f = h5py.File(C.H5, "r")
ct = cat(f, "Cell Type"); st = cat(f, "SampleType"); sid = cat(f, "SampleID")
rk = f["raw"]["var"].attrs.get("_index", b"_index"); rk = rk.decode() if isinstance(rk, bytes) else rk
sym = list(dec(f["raw"]["var"][rk][:]))
syms = [s for s in GENES if s in sym]
gi = [sym.index(s) for s in syms]
rx = f["raw"]["X"]
X = sp.csr_matrix((rx["data"][:], rx["indices"][:], rx["indptr"][:]), shape=tuple(rx.attrs["shape"]))
f.close()
pos = np.asarray(X[:, gi].todense()) > 0

leuk = np.isin(ct, C.LEUK_CELLTYPES) & (st == "AML")
marrow = np.isin(ct, ["0_HSPC", "Myeloid_Pro"])
aml_samples = sorted(set(sid[leuk]))
jc = syms.index("CD96")

cd96_frac = {s: 100 * pos[leuk & (sid == s)][:, jc].mean() for s in aml_samples if (leuk & (sid == s)).sum()}
tgt_pts = [s for s, fr in cd96_frac.items() if fr >= CD96POS]

def metrics(extra_idx):
    cov = pos[:, jc] if extra_idx is None else (pos[:, jc] | pos[:, extra_idx])
    within = np.mean([100 * cov[leuk & (sid == s)].mean() for s in tgt_pts])
    perpt = [100 * cov[leuk & (sid == s)].mean() for s in aml_samples if (leuk & (sid == s)).sum()]
    between = 100 * np.mean([p >= COVER for p in perpt])
    tox = 100 * cov[marrow].mean()
    return within, between, tox

rows = []
w0, b0, t0 = metrics(None)
rows.append(dict(partner="CD96 alone", label="CD96 alone", group="lead", within=round(w0, 2),
                 between=round(b0, 2), tox=round(t0, 2), highly_toxic=False, base=True))
for g in syms:
    if g == "CD96": continue
    disp = GENES[g]; w, b, t = metrics(syms.index(g))
    rows.append(dict(partner=g, label="+" + disp, group=GROUP[disp], within=round(w, 2),
                     between=round(b, 2), tox=round(t, 2), highly_toxic=disp in HIGHLY_TOXIC, base=False))
df = pd.DataFrame(rows).sort_values("within", ascending=False)
for d in (SD, SD08):
    df.to_csv(os.path.join(d, "fig7_combo.csv"), index=False)
print("CD96-targetable patients (>=10%% CD96+): %d/%d | corrected leukemic=%s & AML\n" % (len(tgt_pts), len(aml_samples), C.LEUK_CELLTYPES))
print(df[["label", "group", "within", "between", "tox", "highly_toxic"]].to_string(index=False))
