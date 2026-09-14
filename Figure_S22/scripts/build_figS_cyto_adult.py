#!/usr/bin/env python
# ====================================================================
# build_figS_cyto_adult.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_figS_cyto_adult.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Adult twin of the cytogenetic-faceted supplement — BeatAML (aml_ohsu_2018) bulk, the same
8 panel-F targets violined, faceted by adult cytogenetic subtype.

Cytogenetic subtype per sample is classified from the BeatAML clinical FUSION call, falling
back to the KARYOTYPE string for the no-fusion majority:
  PML-RARA (APL) | CBF [t(8;21)+inv(16)] | KMT2Ar | Normal karyotype | Other
Joins to the panel-F per-sample log2(CPM+1) values (fig5F_adult_beataml_violin.csv).
Facets restricted to subtypes with >=MIN_N samples.
-> source_data/figS_cyto_adult_violin.csv  (sample, gene, log2cpm, label, group, cyto)
"""
import os, re, json, urllib.request, collections
import pandas as pd

API = "https://www.cbioportal.org/api"
HERE = os.path.dirname(os.path.abspath(__file__))
SD = os.path.join(HERE, "..", "source_data")
STUDY = "aml_ohsu_2018"
MIN_N = 5   # include CBF/KMT2Ar (n=5); render shows small facets as points+median, not violins
NORMAL_RE = re.compile(r"^4[5-7],X[XY](\[\d+\])?$")   # single normal clone (46,XX/46,XY[n])


def _get(u):
    return json.load(urllib.request.urlopen(urllib.request.Request(
        u, headers={"Accept": "application/json"}), timeout=90))


def classify(fusion, karyo):
    f = (fusion or "").strip()
    if f == "PML-RARA":
        return "PML-RARA"
    if f in ("CBFB-MYH11", "RUNX1-RUNX1T1"):
        return "CBF [t(8;21)/inv(16)]"
    if "KMT2A" in f:
        return "KMT2Ar"
    if f and f not in ("None", "Unknown"):
        return "Other"
    # no/unknown fusion -> use karyotype
    k = (karyo or "").strip()
    if k and NORMAL_RE.match(k.replace(" ", "")):
        return "Normal karyotype"
    return "Other"


vio = pd.read_csv(os.path.join(SD, "fig5F_adult_beataml_violin.csv")).rename(columns={"log2": "log2cpm"})
# include ALL panel-F genes (now 12: + 4 PPAC-specific) so the supplement tracks the main panel
cohort = set(vio["sample"].unique())

cl = _get(f"{API}/studies/{STUDY}/clinical-data?clinicalDataType=SAMPLE&pageSize=200000")
attr = collections.defaultdict(dict)
for r in cl:
    if r["sampleId"] in cohort and r["clinicalAttributeId"] in ("FUSION", "KARYOTYPE"):
        attr[r["sampleId"]][r["clinicalAttributeId"]] = r["value"]
s2c = {s: classify(a.get("FUSION"), a.get("KARYOTYPE")) for s, a in attr.items()}
vio["cyto"] = vio["sample"].map(s2c)
vio = vio.dropna(subset=["cyto"])

n_by_cyto = vio.groupby("cyto")["sample"].nunique().sort_values(ascending=False)
keep = n_by_cyto[n_by_cyto >= MIN_N].index.tolist()
out = vio[vio["cyto"].isin(keep)].copy()
out.to_csv(os.path.join(SD, "figS_cyto_adult_violin.csv"), index=False)

print("adult BeatAML cohort, cytogenetic subtype sample counts:")
print(n_by_cyto.to_string())
print(f"\n  kept (>= {MIN_N}): {keep}")
print(f"  dropped (< {MIN_N}): {dict(n_by_cyto[n_by_cyto < MIN_N])}")
print(f"  wrote figS_cyto_adult_violin.csv  ({out['sample'].nunique()} samples x "
      f"{out['gene'].nunique()} genes x {len(keep)} subtypes)")
