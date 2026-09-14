#!/usr/bin/env python
"""Figure 6 panels D (radar profiles) + E (cumulative coverage), corrected scheme. Recomputed from the
MAIN atlas raw counts (h5py) on the corrected leukemic definition (Cell_Type in config.LEUK_CELLTYPES
AND SampleType=="AML"), normal HSPC/Myeloid_Pro + mature-immune from Cell_Type, so it is consistent
with the rest of Figures 5-7.
D: per-target metric profile (all axes 0-100, higher = better): AML coverage, patient breadth, marrow
   sparing, mature-immune sparing, vital-organ sparing (worst normal-tissue cell type from CZ Census).
E: CD96-anchored greedy cumulative % of AML patients covered (>=20% of their blasts positive).
-> source_data/fig6D_radar.csv, fig5E_coverage.csv (+ 08_Source_Data mirror). No fabricated values."""
import os, sys, numpy as np, pandas as pd, scipy.sparse as sp, h5py
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C

SD = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/shared_CD96_pipeline/source_data"
SD08 = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_5"
SEL = {"CD96": "CD96", "CD9": "CD9", "SUCNR1": "SUCNR1", "IL1RAP": "IL1RAP", "FLT3": "FLT3",
       "CD33": "CD33", "IL3RA": "CD123", "CLEC12A": "CLL-1", "CD7": "CD7", "TNFRSF4": "TNFRSF4",
       "ABCA7": "ABCA7", "ITGAX": "ITGAX"}                        # symbol -> display
GROUP = {"CD96": "lead", "CD9": "lead", "SUCNR1": "lead", "IL1RAP": "lead", "FLT3": "control",
         "CD33": "clinical", "IL3RA": "clinical", "CLEC12A": "clinical", "CD7": "ppac",
         "TNFRSF4": "ppac", "ABCA7": "ppac", "ITGAX": "ppac"}
SELECTED = ["CD96", "SUCNR1", "CD33", "CLEC12A", "CD7", "TNFRSF4", "ITGAX"]   # E coverage
COVER_THRESH = 20.0
MATURE = ["CD14_Monocyte", "CD16_Monocyte", "mDC", "pDC", "Macrophage", "Naïve_CD4T",
          "Naïve_CD8T", "CTL", "Activated_CD4T", "NK", "CD20+B", "PlasmaB"]
VITAL_KW = ("cardiac muscle", "cardiomyocyte", "neuron", "hepatocyte", "kidney", "renal",
            "nephron", "podocyte", "pneumocyte", "alveolar", "endothelial cell")

dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
def cat(f, k):
    o = f["obs"][k]; return dec(o["categories"][:])[o["codes"][:].astype(int)]

f = h5py.File(C.H5, "r")
ct = cat(f, "Cell Type"); st = cat(f, "SampleType"); sid = cat(f, "SampleID")
rk = f["raw"]["var"].attrs.get("_index", b"_index"); rk = rk.decode() if isinstance(rk, bytes) else rk
sym = list(dec(f["raw"]["var"][rk][:]))
syms = [s for s in SEL if s in sym]
gi = [sym.index(s) for s in syms]
rx = f["raw"]["X"]
X = sp.csr_matrix((rx["data"][:], rx["indices"][:], rx["indptr"][:]), shape=tuple(rx.attrs["shape"]))
f.close()
pos = np.asarray(X[:, gi].todense()) > 0                         # cells x 12 genes

leuk = np.isin(ct, C.LEUK_CELLTYPES) & (st == "AML")            # corrected leukemic set
hspc = ct == "0_HSPC"; mye = ct == "Myeloid_Pro"; mature = np.isin(ct, MATURE)
aml_samples = sorted(set(sid[leuk]))

patmat = np.zeros((len(aml_samples), len(syms)))
for i, s in enumerate(aml_samples):
    m = leuk & (sid == s)
    if m.sum(): patmat[i] = 100.0 * pos[m].mean(0)

cc = pd.read_csv(os.path.join(SD, "fig6tox_celltypes.csv"))     # CZ Census (independent, unchanged)
lc = cc["cell_type"].str.lower()
vital = cc[lc.apply(lambda n: any(k in n for k in VITAL_KW)) & (cc["n"] >= 1000)]
organ_worst = {SEL[s]: float(vital.loc[vital.label == SEL[s], "pct"].max() or 0) for s in syms}

rows = []
for s in syms:
    j = syms.index(s)
    aml = 100.0 * pos[leuk][:, j].mean()
    breadth = 100.0 * (patmat[:, j] >= COVER_THRESH).mean()
    marrow_leak = max(100.0 * pos[hspc][:, j].mean(), 100.0 * pos[mye][:, j].mean())
    mature_leak = 100.0 * pos[mature][:, j].mean()
    rows.append(dict(gene=s, label=SEL[s], group=GROUP[s], AML_coverage=round(aml, 1),
                     Patient_breadth=round(breadth, 1), Marrow_sparing=round(100 - marrow_leak, 1),
                     Immune_sparing=round(100 - mature_leak, 1),
                     Organ_sparing=round(100 - organ_worst[SEL[s]], 1)))
radar = pd.DataFrame(rows)

covered = (patmat >= COVER_THRESH)
cd96j = syms.index("CD96"); got = covered[:, cd96j].copy(); chosen = [cd96j]
remaining = set(syms.index(s) for s in SELECTED) - {cd96j}
steps = [dict(step=1, gene="CD96", label="CD96", group="lead", added=int(covered[:, cd96j].sum()),
              cum_pct=round(100.0 * got.mean(), 1), n_patients=len(aml_samples))]
while remaining:
    best, best_new = None, -1
    for j in remaining:
        new = int((covered[:, j] & ~got).sum())
        if new > best_new: best_new, best = new, j
    chosen.append(best); remaining.discard(best); got |= covered[:, best]
    steps.append(dict(step=len(chosen), gene=syms[best], label=SEL[syms[best]], group=GROUP[syms[best]],
                      added=best_new, cum_pct=round(100.0 * got.mean(), 1), n_patients=len(aml_samples)))

for d in (SD, SD08):
    radar.to_csv(os.path.join(d, "fig6D_radar.csv"), index=False)
    pd.DataFrame(steps).to_csv(os.path.join(d, "fig5E_coverage.csv"), index=False)
print("AML patients: %d | corrected leukemic=%s & AML | threshold >=%.0f%%\n" % (len(aml_samples), C.LEUK_CELLTYPES, COVER_THRESH))
print("=== D radar ===\n", radar.to_string(index=False))
print("\n=== E greedy coverage ===\n", pd.DataFrame(steps)[["step", "label", "added", "cum_pct"]].to_string(index=False))
