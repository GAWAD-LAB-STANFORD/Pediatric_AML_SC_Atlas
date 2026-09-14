#!/usr/bin/env python
"""Every 2-target combination of the 12 antigens (66 pairs) + 12 singles — toxicity vs efficacy, within
AND across patients. Corrected scheme: recomputed from the MAIN atlas raw counts (h5py) on the corrected
leukemic definition (Cell_Type in config.LEUK_CELLTYPES AND SampleType=="AML"), normal marrow =
Cell_Type 0_HSPC/Myeloid_Pro. For each union "X OR Y" per cell:
  toxicity = % of normal marrow cells positive; within = % blasts covered (mean across AML patients);
  across   = % of AML patients with >=20% of blasts covered.
-> source_data/fig7_allpairs.csv (+ 08_Source_Data mirror). No fabricated values."""
import os, sys, itertools, numpy as np, pandas as pd, scipy.sparse as sp, h5py
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C

SD = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/shared_CD96_pipeline/source_data"
SD08 = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_5"
GENES = {"CD96": "CD96", "CD9": "CD9", "SUCNR1": "SUCNR1", "IL1RAP": "IL1RAP", "FLT3": "FLT3",
         "CD33": "CD33", "IL3RA": "CD123", "CLEC12A": "CLL-1", "CD7": "CD7", "TNFRSF4": "TNFRSF4",
         "ABCA7": "ABCA7", "ITGAX": "ITGAX"}
COVER = 20.0

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
aml = sorted(set(sid[leuk]))
leuk_s = {s: (leuk & (sid == s)) for s in aml if (leuk & (sid == s)).sum()}

def metrics(cols):
    cov = pos[:, cols[0]].copy()
    for c in cols[1:]: cov |= pos[:, c]
    perpt = [100 * cov[m].mean() for m in leuk_s.values()]
    return (round(100 * cov[marrow].mean(), 2), round(float(np.mean(perpt)), 2),
            round(100 * float(np.mean([p >= COVER for p in perpt])), 2))

rows = []
for g in syms:
    tox, wi, ac = metrics([syms.index(g)])
    rows.append(dict(combo=GENES[g], g1=GENES[g], g2="", kind="single",
                     cd96=int(g == "CD96"), toxicity=tox, within=wi, across=ac))
for a_, b_ in itertools.combinations(syms, 2):
    da, db = GENES[a_], GENES[b_]
    tox, wi, ac = metrics([syms.index(a_), syms.index(b_)])
    rows.append(dict(combo=f"{da}+{db}", g1=da, g2=db, kind="pair",
                     cd96=int("CD96" in (da, db)), toxicity=tox, within=wi, across=ac))
df = pd.DataFrame(rows)
for d in (SD, SD08):
    df.to_csv(os.path.join(d, "fig7_allpairs.csv"), index=False)
print("wrote fig7_allpairs.csv (%d rows: %d singles + 66 pairs) | corrected leukemic=%s & AML\n"
      % (len(df), len(syms), C.LEUK_CELLTYPES))
print(df[df.kind == "pair"].sort_values("within", ascending=False).head(10)
      [["combo", "cd96", "toxicity", "within", "across"]].to_string(index=False))
