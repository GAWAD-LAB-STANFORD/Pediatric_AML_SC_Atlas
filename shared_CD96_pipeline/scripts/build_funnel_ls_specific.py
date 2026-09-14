#!/usr/bin/env python
# ====================================================================
# build_funnel_ls_specific.py  |  Upstream builder (regenerates source_data)
# Poor-prognosis-STATE (LS) surface-target discovery -- replaces the retired
# per-PPAC version (build_funnel_ppac_specific.py) on the corrected 49-state scheme.
# ====================================================================
"""Poor-prognosis leukemic-STATE (LS) surface-target discovery.

Ranks surface markers WITHIN each poor-prognosis leukemic state (LS) by the % of
that state's cells positive (count > 0), behind the same toxicity funnel as the
retired PPAC version: normal HSPC <= GATE% AND normal Myeloid_Pro <= GATE%, and the
TCR/CD3 family dropped (MPAL-like co-expression). Uses the corrected 49-state
assignment (reference_assignments res4.50 -> config.stable_raws; the "49 file").

Groups (parallel to the old pooled + PPAC_1..5):
  "All poor LS"  = all 8 poor-prognosis LS pooled
  the 5 MULTI-patient poor LS (generalizable): LS_10, LS_19, LS_20, LS_21, LS_27
The 3 single-patient poor LS (LS_23/LS_38/LS_41; >=77% of cells from one patient)
are pooled into "All poor LS" but NOT given individual panels (single-patient
artifacts, not generalizable) -- an improvement over the PPAC version, which had to
show a single-patient PPAC_5. Per-group patient composition is audited so any
single-patient signal is flagged. No fabricated values.

Outputs:
  source_data/funnel_ls_specific.csv               (long: group x gated-gene metrics)
  source_data/funnel_ls_specific_composition.csv   (per-group patient composition)
"""
import os, sys
import numpy as np, pandas as pd
import h5py
from scipy.sparse import csr_matrix

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3")
import config as C
SURF = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/ref/surfaceome_genes.txt"
PROG = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3/LS_wholecohort_prognosis.csv"
OUT     = os.path.join(HERE, "..", "source_data", "funnel_ls_specific.csv")
OUT_COMP= os.path.join(HERE, "..", "source_data", "funnel_ls_specific_composition.csv")

GATE_HSPC, GATE_MYE = 10.0, 10.0      # toxicity tolerance (matches the PPAC version)
REC_THRESH = 0.20
MIN_CELLS_PER_PT = 20
TOP_N = 12
SINGLE_PT_FRAC = 0.77                  # >= this share from one patient => single-patient state
TCR_CD3_EXCLUDE = {"CD3D","CD3E","CD3G","CD247","TRAC","TRBC1","TRBC2","TRGC1","TRGC2","TRDC"}
CLIN = {"CD96","CD33","IL3RA","CLEC12A","FLT3","MSLN","CD70","ADGRG1","IL1RAP"}
GLABEL = {"IL3RA":"CD123","CLEC12A":"CLL-1","ADGRG1":"GPR56"}

# ---- poor-prognosis LS (from the 49-state whole-cohort prognosis) ----
prog = pd.read_csv(PROG)
poor = prog[prog.prognosis == "poor"].copy()
poor["single_pt"] = poor.top_pt_frac >= SINGLE_PT_FRAC
poor_all   = poor.LS.tolist()
poor_multi = poor[~poor.single_pt].sort_values("n_cells", ascending=False).LS.tolist()
print(f"poor LS (n={len(poor_all)}): {poor_all}")
print(f"  multi-patient (panels): {poor_multi}")
print(f"  single-patient (pooled only): {poor[poor.single_pt].LS.tolist()}")

# ---- load atlas: barcodes, Cell Type (gate), SampleID (patients), raw counts ----
dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
f = h5py.File(C.H5, "r")
obs = dec(f["obs"]["_index"][:])
def catcol(name):
    g = f["obs"][name]
    return pd.Series(dec(g["categories"][:])[g["codes"][:]], index=obs) if isinstance(g, h5py.Group) \
           else pd.Series(dec(g[:]), index=obs)
ct   = catcol("Cell Type")
samp = catcol("SampleID")
rv   = dec(f["raw"]["var"]["_index"][:]); rpos = {g: i for i, g in enumerate(rv)}
surf = [g.strip() for g in open(SURF) if g.strip()]
present = np.array([g for g in surf if g in rpos])
gi = [rpos[g] for g in present]
Xr = f["raw"]["X"]; shape = tuple(Xr.attrs["shape"]) if "shape" in Xr.attrs else (len(obs), len(rv))
X = csr_matrix((Xr["data"][:], Xr["indices"][:], Xr["indptr"][:]), shape=shape)
pos = (np.asarray(X[:, gi].todense()) > 0)          # cells x Npresent
f.close()
print(f"{len(present)} surfaceome genes detectable in scRNA raw counts")

# ---- LS assignment (the corrected "49 file") ----
ra = pd.read_csv(C.REFA); ra["barcode"] = ra["barcode"].astype(str).str.replace(r"^b'|'$", "", regex=True)
raws, ls = C.stable_raws(); raw2ls = {rw: ls[i] for i, rw in enumerate(raws)}
bc2ls = dict(zip(ra.barcode, ra["res4.50"].map(raw2ls)))
lscol = pd.Series([bc2ls.get(b) for b in obs], index=obs)

# ---- toxicity funnel (HSPC<=10%, Myeloid<=10%, drop TCR/CD3) ----
hspc = (ct == "0_HSPC").values
mye  = (ct == "Myeloid_Pro").values
hspc_pct = 100 * pos[hspc].mean(0)
mye_pct  = 100 * pos[mye].mean(0)
is_tcr = np.array([g in TCR_CD3_EXCLUDE for g in present])
gate = (hspc_pct <= GATE_HSPC) & (mye_pct <= GATE_MYE) & (~is_tcr)
print(f"eligible genes after funnel: {gate.sum()}")

groups = [("All poor LS", lscol.isin(poor_all).values)] + \
         [(g, (lscol == g).values) for g in poor_multi]

comp_rows, long_rows = [], []
for gname, gmask in groups:
    n_cells = int(gmask.sum())
    pt_counts = samp[gmask].value_counts()
    pt_counts = pt_counts[pt_counts > 0]
    contrib = pt_counts[pt_counts >= MIN_CELLS_PER_PT].index.tolist()
    top_share = 100 * pt_counts.iloc[0] / n_cells if n_cells else 0.0
    comp_rows.append(dict(group=gname, n_cells=n_cells, n_patients_any=int((pt_counts > 0).sum()),
                          n_patients_contrib=len(contrib),
                          dominant_patient_share_pct=round(top_share, 1),
                          single_patient_flag=(len(contrib) < 3)))
    cov = 100 * pos[gmask].mean(0)
    if contrib:
        hit = np.zeros(len(present), int)
        for s in contrib:
            m = gmask & (samp == s).values
            hit += (pos[m].mean(0) >= REC_THRESH)
        rec = 100 * hit / len(contrib)
    else:
        rec = np.zeros(len(present))
    for i in np.where(gate)[0]:
        g = present[i]
        long_rows.append(dict(group=gname, gene=g, label=GLABEL.get(g, g),
                              cov_pct=round(float(cov[i]), 2), rec_pct=round(float(rec[i]), 2),
                              HSPC=round(float(hspc_pct[i]), 2), Myeloid=round(float(mye_pct[i]), 2),
                              n_contrib_patients=len(contrib), clinical=g in CLIN))

comp = pd.DataFrame(comp_rows); long = pd.DataFrame(long_rows)
long.to_csv(OUT, index=False); comp.to_csv(OUT_COMP, index=False)
pd.set_option("display.width", 160)
print("\n=== PER-LS PATIENT COMPOSITION ===")
print(comp.to_string(index=False))
for gname, _ in groups:
    sub = long[long.group == gname].sort_values("cov_pct", ascending=False).head(TOP_N)
    print(f"\n=== {gname} ===")
    print(sub[["label","cov_pct","rec_pct","HSPC","Myeloid","clinical"]].to_string(index=False))
