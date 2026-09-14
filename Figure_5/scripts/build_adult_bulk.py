#!/usr/bin/env python
# ====================================================================
# build_adult_bulk.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_adult_bulk.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Pull adult-AML bulk expression for the panel-E genes from cBioPortal (BeatAML + TCGA-LAML),
restrict to >=80% BM-blast samples, and write violin-ready tables that mirror panel E.
  BeatAML (aml_ohsu_2018): CPM profile -> log2(CPM+1), blast = BM_BLAST_PERCENTAGE
  TCGA-LAML (laml_tcga):   RSEM profile -> log2(RSEM+1), blast = BLAST_COUNT
Outputs: source_data/fig6F_adult_{beataml,tcga}_violin.csv (+ _meta.csv).
"""
import os, json, urllib.request
import numpy as np, pandas as pd

API = "https://www.cbioportal.org/api"
HERE = os.path.dirname(os.path.abspath(__file__))
SD = os.path.join(HERE, "..", "source_data")          # feeds the figure deck (BeatAML = panel F)
ADULT = os.path.join(HERE, "..", "external_data", "adult_aml")   # archive (TCGA: n=29 @>=80%, RSEM)
os.makedirs(ADULT, exist_ok=True)
BLAST_MIN = 80.0
TARGETABLE = 5.0
# symbol -> (entrez, display, group)
GENES = {"CD96": (10225, "CD96", "lead"), "CD9": (928, "CD9", "lead"),
         "SUCNR1": (56670, "SUCNR1", "lead"), "IL1RAP": (3556, "IL1RAP", "lead"),
         "FLT3": (2322, "FLT3", "control"), "CD33": (945, "CD33", "clinical"),
         "IL3RA": (3563, "CD123", "clinical"), "CLEC12A": (160364, "CLL-1", "clinical"),
         # PPAC-specific markers (S17 group; purple) — mark individual PPAC subsets, not pan-AML
         "CD7": (924, "CD7", "ppac"), "TNFRSF4": (7293, "TNFRSF4", "ppac"),
         "ABCA7": (10347, "ABCA7", "ppac"), "ITGAX": (3687, "ITGAX", "ppac")}
ENTREZ2SYM = {v[0]: k for k, v in GENES.items()}
# scRNA normal-marrow toxicity (max %HSPC, %Myeloid_Pro positive) — cohort-independent
# property of the gene vs our healthy-BM scRNA; carried so panel F reads as a twin of E.
HSPC_MYE = {"CD96": 3.02, "CD9": 7.61, "SUCNR1": 7.12, "IL1RAP": 9.65,
            "FLT3": 53.06, "CD33": 28.71, "IL3RA": 12.62, "CLEC12A": 29.67,
            "CD7": 9.09, "TNFRSF4": 0.56, "ABCA7": 9.77, "ITGAX": 6.63}

def _get(url):
    return json.load(urllib.request.urlopen(urllib.request.Request(
        url, headers={"Accept": "application/json"}), timeout=90))
def _post(url, body):
    return json.load(urllib.request.urlopen(urllib.request.Request(
        url, data=json.dumps(body).encode(), method="POST",
        headers={"Content-Type": "application/json", "Accept": "application/json"}), timeout=180))

def blast_by_sample(study, attr):
    """map sampleId -> blast value, resolving patient-level attrs via patientId."""
    samp = _get(f"{API}/studies/{study}/clinical-data?clinicalDataType=SAMPLE&pageSize=100000")
    pat = _get(f"{API}/studies/{study}/clinical-data?clinicalDataType=PATIENT&pageSize=100000")
    s_blast = {d["sampleId"]: d["value"] for d in samp if d["clinicalAttributeId"] == attr}
    p_blast = {d["patientId"]: d["value"] for d in pat if d["clinicalAttributeId"] == attr}
    s2p = {d["sampleId"]: d["patientId"] for d in samp}
    out = {}
    for s in s2p:
        v = s_blast.get(s) or p_blast.get(s2p[s])
        try:
            out[s] = float(v)
        except (TypeError, ValueError):
            pass
    return out

def fetch(label, study, profile, blast_attr, logbase, outdir=SD):
    print(f"\n=== {label} ({study}) ===")
    blast = blast_by_sample(study, blast_attr)
    keep = {s for s, b in blast.items() if b >= BLAST_MIN}
    print(f"  samples with blast metadata: {len(blast)};  >= {BLAST_MIN:.0f}% blasts: {len(keep)}")
    mol = _post(f"{API}/molecular-profiles/{profile}/molecular-data/fetch",
                {"entrezGeneIds": [v[0] for v in GENES.values()], "sampleListId": f"{study}_all"})
    df = pd.DataFrame([(m["sampleId"], ENTREZ2SYM.get(m["entrezGeneId"]), m["value"]) for m in mol],
                      columns=["sample", "sym", "value"]).dropna()
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df[df["sample"].isin(keep)].dropna()
    df["log2"] = np.log2(df["value"].clip(lower=0) + 1.0)
    df["label"] = df["sym"].map(lambda s: GENES[s][1]); df["group"] = df["sym"].map(lambda s: GENES[s][2])
    rows, meta = [], []
    for sym, (ent, lab, grp) in GENES.items():
        v = df.loc[df.sym == sym, "log2"]
        if not len(v):
            continue
        for s, val in zip(df.loc[df.sym == sym, "sample"], v):
            rows.append(dict(sample=s, gene=sym, log2=float(val), label=lab, group=grp))
        meta.append(dict(gene=sym, label=lab, group=grp, n=int(len(v)),
                         median=round(float(v.median()), 3),
                         pct_targetable=round(100 * float((v > TARGETABLE).mean()), 1),
                         hspc_mye=HSPC_MYE[sym], toxic=bool(HSPC_MYE[sym] > 10)))
    suf = label.lower()
    pd.DataFrame(rows).to_csv(os.path.join(outdir, f"fig6F_adult_{suf}_violin.csv"), index=False)
    m = pd.DataFrame(meta).sort_values("median", ascending=False)
    m.to_csv(os.path.join(outdir, f"fig6F_adult_{suf}_meta.csv"), index=False)
    print(f"  n samples used: {df['sample'].nunique()};  expression units: {logbase}")
    print(m.to_string(index=False))
    return m

# BeatAML = the clean adult twin (CPM, n=104) -> source_data -> panel F.
fetch("BeatAML", "aml_ohsu_2018", "aml_ohsu_2018_mrna_seq_cpm", "BM_BLAST_PERCENTAGE", "log2(CPM+1)")
# TCGA-LAML = downloaded for completeness but archived, NOT in the deck: only 29 samples reach
# >=80% blasts and cBioPortal serves RSEM (different scale, not CPM-comparable to BeatAML/TARGET).
fetch("TCGA", "laml_tcga", "laml_tcga_rna_seq_v2_mrna", "BLAST_COUNT", "log2(RSEM+1)", outdir=ADULT)
