#!/usr/bin/env python
# ====================================================================
# build_figS_cyto_pediatric.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_figS_cyto_pediatric.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Supplement to Figure-5 panel E — the 8 panel-E targets violined, FACETED by cytogenetic
subtype, in the PEDIATRIC TARGET bulk cohort (>=80% blasts, n=475).

Joins the panel-E per-sample log2(CPM+1) values (source_data/fig5F_violin.csv) to each
sample's cytogenetic subtype, classified from the TARGET clinical table
(Figure7F_TARGET_CD96_merged.csv) with the same fusion/karyotype rule used in the CD96 work.
Facets are restricted to subtypes with >=MIN_N samples (violins of n<8 are not shown).
-> source_data/figS_cyto_pediatric_violin.csv  (sample, gene, log2cpm, label, group, cyto)
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SD = os.path.join(ROOT, "source_data")
MERGED = os.path.abspath(os.path.join(ROOT, "..", "Figure7", "Figure7F_TARGET_CD96_merged.csv"))
MIN_N = 8


def fam(row):
    """TARGET sample -> cytogenetic subtype (fusion first, else primary cytogenetic code)."""
    fz = str(row.get("Gene Fusion", "")).strip()
    if fz and fz.lower() != "nan":
        if fz.startswith("KMT2A-") or fz == "other MLL":
            return "KMT2Ar"
        if fz == "RUNX1-RUNX1T1":
            return "t(8;21)"
        if fz == "CBFB-MYH11":
            return "inv(16)"
        if fz.startswith("NUP98-"):
            return "NUP98r"
        if fz == "CBFA2T3-GLIS2":
            return "CBFA2T3::GLIS2"
        if fz == "MLLT10-PICALM":
            return "PICALM::MLLT10"
        if fz == "DEK-NUP214":
            return "DEK::NUP214"
        if fz == "NPM1-MLF1":
            return "NPM1::MLF1"
        return "Other"
    return {"Normal": "CN", "t(8;21)": "t(8;21)", "inv(16)": "inv(16)",
            "PML-RARA": "PML-RARA", "MLL": "KMT2Ar"}.get(
        str(row.get("Primary Cytogenetic Code", "")).strip(), "Other")


vio = pd.read_csv(os.path.join(SD, "fig5F_violin.csv"))         # panel-E genes x 475 samples, log2cpm
# include ALL panel-E genes (4 leads + FLT3 + 3 clinical + 4 PPAC-specific) so the supplement
# tracks the main panel exactly
m = pd.read_csv(MERGED)
m["cyto"] = m.apply(fam, axis=1)
s2c = dict(zip(m["sample"], m["cyto"]))
vio["cyto"] = vio["sample"].map(s2c)

miss = vio["cyto"].isna().sum()
if miss:
    print(f"  WARNING: {miss} rows had no cytogenetic match (left as NA)")
vio = vio.dropna(subset=["cyto"])

# subtype sample counts (per subtype, not per row)
n_by_cyto = vio.groupby("cyto")["sample"].nunique().sort_values(ascending=False)
keep = n_by_cyto[n_by_cyto >= MIN_N].index.tolist()
dropped = n_by_cyto[n_by_cyto < MIN_N]
out = vio[vio["cyto"].isin(keep)].copy()
out.to_csv(os.path.join(SD, "figS_cyto_pediatric_violin.csv"), index=False)

print("pediatric TARGET cohort, subtype sample counts:")
print(n_by_cyto.to_string())
print(f"\n  kept (>= {MIN_N}): {keep}")
print(f"  dropped (< {MIN_N}, too few for a violin): {dict(dropped)}")
print(f"  wrote figS_cyto_pediatric_violin.csv  ({out['sample'].nunique()} samples x "
      f"{out['gene'].nunique()} genes x {len(keep)} subtypes)")
