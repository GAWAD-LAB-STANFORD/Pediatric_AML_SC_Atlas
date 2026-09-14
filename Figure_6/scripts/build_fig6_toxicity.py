#!/usr/bin/env python
# ====================================================================
# build_fig6_toxicity.py  |  Figure 6 (builder)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_fig6_toxicity.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""NEW Figure 6 (toxicity disqualification) — body-wide normal expression of the 12 panel-5E
genes from the fresh CZ Census pull (external_data/census/).

  A organ atlas       -> source_data/fig6tox_organ.csv
        per organ, aggregated over PARENCHYMAL (non-haematopoietic) cell types only, so a
        lymphoid gene (e.g. CD7) does not light up every organ via infiltrating immune cells.
  B haematopoiesis    -> source_data/fig6tox_hema.csv
        census-wide normal haematopoietic cell types (HSPC -> myeloid/lymphoid/erythroid/mega),
        the on-target marrow-toxicity axis for an AML target.

% positive = fraction of normal 10x cells with >=1 read; mean = mean log1p(CP10k). Per-gene
colour scaling is done in R so no single broad gene (CD9) dominates the shared scale.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
CEN = os.environ.get("CENSUS_DIR", os.path.join(ROOT, "external_data", "census"))
SD = os.path.join(ROOT, "source_data")

# gene -> (display label, group).  group: lead | control | clinical | ppac
GENES = {"CD96": ("CD96", "lead"), "CD9": ("CD9", "lead"), "SUCNR1": ("SUCNR1", "lead"),
         "IL1RAP": ("IL1RAP", "lead"), "FLT3": ("FLT3", "control"), "CD33": ("CD33", "clinical"),
         "IL3RA": ("CD123", "clinical"), "CLEC12A": ("CLL-1", "clinical"),
         "CD7": ("CD7", "ppac"), "TNFRSF4": ("TNFRSF4", "ppac"), "ABCA7": ("ABCA7", "ppac"),
         "ITGAX": ("ITGAX", "ppac")}
SYMS = list(GENES)

# haematopoietic / immune cell-type keywords (lowercase substring match on the CL term).
HEMA_KW = ("t cell", "t-cell", "t-helper", "t helper", "cd4-positive", "cd8-positive", "cd4 ", "cd8 ",
           "regulatory t", "helper t", "follicular helper", "cytotoxic t", "effector t", "memory t",
           "naive t", "alpha-beta t", "double-positive", "double-negative", "th1", "th2", "th17",
           "gamma-delta", "mait", "mucosal invariant", "thymocyte", "b cell", "pro-b", "pre-b",
           "plasma cell", "plasmablast", "memory b", "naive b", "germinal", "natural killer", "nk cell",
           "innate lymphoid", "lymphocyte", "leukocyte", "monocyte", "macrophage", "kupffer",
           "microglial", "dendritic", "langerhans", "myeloid cell", "myeloid leukocyte", "hofbauer",
           "granulocyte", "neutrophil", "eosinophil", "basophil", "mast cell",
           "erythroid", "erythrocyte", "normoblast", "reticulocyte", "erythroblast",
           "megakaryocyte", "platelet", "haematopoietic", "hematopoietic", "myeloblast",
           "promyelocyte", "myelocyte", "lymphoid", "plasmacytoid", "phagocyte", "mononuclear",
           "histiocyte", "kupffer", "osteoclast", "blood cell", "antigen presenting",
           "antigen-presenting", "myeloid lineage", "common myeloid", "common lymphoid",
           "t-regulatory", "myeloid progenitor", "lymphoid progenitor")
# ordered haematopoietic hierarchy for panel B (substring -> canonical bucket, in order)
HEMA_ORDER = [
    ("hematopoietic stem", "HSC"), ("haematopoietic stem", "HSC"),
    ("common myeloid", "CMP"), ("granulocyte monocyte progenitor", "GMP"),
    ("common lymphoid", "CLP"), ("megakaryocyte-erythroid", "MEP"),
    ("myeloblast", "myeloblast"), ("promyelocyte", "promyelocyte"),
    ("neutrophil", "neutrophil"), ("granulocyte", "granulocyte"),
    ("classical monocyte", "monocyte"), ("monocyte", "monocyte"),
    ("macrophage", "macrophage"), ("conventional dendritic", "cDC"),
    ("plasmacytoid dendritic", "pDC"), ("dendritic", "DC"),
    ("erythroid", "erythroid"), ("erythrocyte", "erythroid"),
    ("megakaryocyte", "megakaryocyte"),
    ("naive b", "naive B"), ("memory b", "memory B"), ("b cell", "B cell"),
    ("plasma cell", "plasma cell"),
    ("cd4-positive", "CD4 T"), ("cd8-positive", "CD8 T"), ("regulatory t", "Treg"),
    ("natural killer", "NK"),
]


def is_hema(name):
    n = str(name).lower()
    return any(k in n for k in HEMA_KW)


def to_long(df, colname):
    rows = []
    for _, r in df.iterrows():
        for s in SYMS:
            pc, me = f"{s}_pct", f"{s}_mean"
            if pc in df.columns and pd.notna(r.get(pc)):
                rows.append(dict(gene=s, label=GENES[s][0], group=GENES[s][1],
                                 **{colname: r[colname]}, pct=float(r[pc]),
                                 mean=float(r[me]) if pd.notna(r.get(me)) else np.nan,
                                 n=int(r["n_cells"])))
    return pd.DataFrame(rows)


# ---------- A: parenchymal organ atlas ----------
ct = pd.read_csv(os.path.join(CEN, "census_by_celltype.csv"))
ct = ct[~ct["cell_type"].map(is_hema)].copy()             # parenchymal only
for s in SYMS:                                             # weighted-mean numerator (mean*n)
    ct[f"{s}_mxn"] = ct[f"{s}_mean"] * ct["n_cells"]
g = ct.groupby("tissue")
organ = g["n_cells"].sum().rename("n_cells").to_frame()
for s in SYMS:                                            # vectorized groupby sums (no .apply)
    organ[f"{s}_pct"] = 100.0 * g[f"{s}_pos"].sum() / organ["n_cells"]
    organ[f"{s}_mean"] = g[f"{s}_mxn"].sum() / organ["n_cells"]
organ = organ.reset_index()
organ = organ[organ["n_cells"] >= 2000]                   # enough parenchymal cells
long_o = to_long(organ, "tissue").rename(columns={"tissue": "organ"})
long_o.to_csv(os.path.join(SD, "fig6tox_organ.csv"), index=False)
print(f"A organ atlas: {organ['tissue'].nunique()} organs (parenchymal) x {len(SYMS)} genes "
      f"-> fig6tox_organ.csv")

# ---------- census-wide per-cell-type roll-up (from by_celltype; partial-safe) ----------
# sum positives & n across tissues per cell type -> equivalent to the end-of-run global roll-up.
cg = pd.read_csv(os.path.join(CEN, "census_by_celltype.csv"))
grows = []
for cname, d in cg.groupby("cell_type"):
    n = int(d["n_cells"].sum())
    row = {"cell_type": cname, "n_cells": n}
    for s in SYMS:
        pos = int(d[f"{s}_pos"].sum())
        row[f"{s}_pos"] = pos
        row[f"{s}_pct"] = 100.0 * pos / n if n else np.nan
        row[f"{s}_mean"] = np.average(d[f"{s}_mean"], weights=d["n_cells"]) if n else np.nan
    grows.append(row)
glob_ct = pd.DataFrame(grows)

# ---------- B: haematopoiesis atlas ----------
gl = glob_ct[glob_ct["cell_type"].map(is_hema) & (glob_ct["n_cells"] >= 500)].copy()
# map each cell type to the FIRST matching hierarchy bucket; keep the largest-n cell type/bucket
def bucket(name):
    n = str(name).lower()
    for kw, b in HEMA_ORDER:
        if kw in n:
            return b
    return None
gl["bucket"] = gl["cell_type"].map(bucket)
gl = gl.dropna(subset=["bucket"])
# for each bucket pick the single largest-n representative cell type (cleanest signal)
rep = gl.sort_values("n_cells", ascending=False).drop_duplicates("bucket")
order = [b for _, b in HEMA_ORDER]
seen, bucket_order = set(), []
for b in order:
    if b not in seen and b in set(rep["bucket"]):
        bucket_order.append(b); seen.add(b)
rep = rep.set_index("bucket").loc[bucket_order].reset_index()
long_h = to_long(rep.rename(columns={"bucket": "cell"}), "cell")
long_h.to_csv(os.path.join(SD, "fig6tox_hema.csv"), index=False)
print(f"B haematopoiesis: {rep['bucket'].nunique()} cell types -> fig6tox_hema.csv")
print("  hierarchy order:", bucket_order)

# ---------- C: deadly-toxicity cell types (vital-organ parenchyma + HSPC) ----------
# normal cell types where on-target killing is lethal; aggregate census-wide per bucket.
# clinically-serious normal cell types where on-target killing is mortal/severe; consolidate the
# census's granular subtypes into recognizable buckets. Expanded from the discovered target-
# specific liabilities (bladder, retina, lens, cornea, epidermis, pancreas). HSPC excluded (panel B).
DEADLY = [
    ("Cardiomyocyte", ["cardiac muscle cell", "cardiomyocyte", "cardiac myocyte"]),
    ("CNS neuron", ["neuron"]),
    ("Retina photoreceptor", ["cone cell", "rod cell", "photoreceptor", "retinal pigment epithelial"]),
    ("Lens fiber", ["lens fiber"]),
    ("Corneal epithelium", ["corneal epithelial"]),
    ("Lung pneumocyte", ["pneumocyte", "alveolar epithelial", "alveolar type"]),
    ("Hepatocyte", ["hepatocyte"]),
    ("Pancreas acinar/islet", ["pancreatic acinar", "pancreatic ductal", "pancreatic a cell",
                               "pancreatic d cell", "type b pancreatic", "islet", "acinar cell"]),
    ("Kidney epithelium", ["kidney epithelial", "renal", "podocyte", "nephron", "kidney proximal",
                           "kidney distal", "kidney loop", "kidney collecting", "kidney connecting"]),
    ("Bladder urothelium", ["urothelial", "umbrella cell"]),
    ("Gut epithelium", ["enterocyte", "epithelial cell of small intestine", "m cell of gut",
                        "epithelial cell of large intestine", "intestinal epithelial", "enteroendocrine"]),
    ("Skin keratinocyte", ["keratinocyte", "spinous", "epidermis"]),
    ("Endothelium", ["endothelial cell"]),
]
drows, kept_deadly = [], []
lc = glob_ct["cell_type"].str.lower()
for bname, kws in DEADLY:
    sub = glob_ct[lc.apply(lambda nm: any(k in nm for k in kws)) & (glob_ct["n_cells"] >= 500)]
    if sub.empty:
        print(f"  [C] {bname}: no cell type >=500 cells (skipped)")
        continue
    kept_deadly.append(bname)
    for s in SYMS:
        # WORST-CASE per gene: the single highest-% cell type in this organ-system bucket
        # (pooling would dilute a peak, e.g. CD96 57% in RPE hidden among low retinal cells).
        j = sub[f"{s}_pct"].astype(float).idxmax()
        r = sub.loc[j]
        drows.append(dict(cell=bname, gene=s, label=GENES[s][0], group=GENES[s][1],
                          pct=round(float(r[f"{s}_pct"]), 3), mean=round(float(r[f"{s}_mean"]), 5),
                          worst_ct=str(r["cell_type"]), n=int(r["n_cells"])))
pd.DataFrame(drows).to_csv(os.path.join(SD, "fig6tox_deadly.csv"), index=False)
print(f"C deadly panel: {len(kept_deadly)} vital cell types -> fig6tox_deadly.csv  ({kept_deadly})")

# ---------- C (data-driven): EVERY non-haematopoietic cell type census-wide, per gene, so the
# render can show each target's OWN top toxicity cell types (not a fixed curated list) ----------
para = glob_ct[(~glob_ct["cell_type"].map(is_hema)) & (glob_ct["n_cells"] >= 500)].copy()
prows = []
for _, r in para.iterrows():
    for s in SYMS:
        prows.append(dict(cell_type=r["cell_type"], gene=s, label=GENES[s][0], group=GENES[s][1],
                          pct=round(float(r[f"{s}_pct"]), 3), mean=round(float(r[f"{s}_mean"]), 5),
                          n=int(r["n_cells"])))
pd.DataFrame(prows).to_csv(os.path.join(SD, "fig6tox_celltypes.csv"), index=False)
print(f"C cell-type pool: {len(para)} non-haematopoietic cell types -> fig6tox_celltypes.csv")
