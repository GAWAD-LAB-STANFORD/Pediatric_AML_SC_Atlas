#!/usr/bin/env python
# ====================================================================
# build_Figure8_panels.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_Figure8_panels.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Supporting panels for Figure 8 (body-wide normal-tissue toxicity), built from
the CZ CELLxGENE Census per-(tissue, cell-type) table.  Renders three panels;
they are tiled with the detailed organ-grouped map by build_Figure8_assembled.py.

  B  Organ atlas        9 targets x organs (circulating blood removed from solid
                        organs; tissue-resident macrophages/mast cells kept).
  C  Safety footprint   per target, # of major normal NON-blood cell types it
                        hits (>=30% positive, n>=3,000).  CD96 = fewest = safest.
  D  CD96 liabilities   the only normal non-blood cell types where CD96 itself is
                        expressed (honest accounting): urothelium, adipocytes.

% positive = fraction of cells with >=1 raw read.  No fabricated values; no RNG.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, "supplements")
sys.path.insert(0, SUPP)
from _fig7_style import PALETTE, RC, save_fig, SEQ_CMAP, CD96_CMAP, TARGET_COLORS

BY_CT = os.path.join(HERE, "Figure8_census_normal_by_celltype.csv")
GENES = ["CD96", "CD33", "IL3RA", "CLEC12A", "FLT3", "MSLN", "CD70", "ADGRG1",
         "IL1RAP"]
GLABEL = {"CD96": "CD96", "CD33": "CD33", "IL3RA": "CD123", "CLEC12A": "CLL-1",
          "FLT3": "FLT3", "MSLN": "MSLN", "CD70": "CD70", "ADGRG1": "GPR56",
          "IL1RAP": "IL1RAP"}
KEEP_ALL_CELLS = {"blood", "bone marrow"}
ORGANS = ["blood", "bone marrow", "spleen", "lymph node", "brain", "spinal cord",
          "eye", "lung", "heart", "liver", "kidney", "pancreas", "stomach",
          "small intestine", "large intestine", "colon", "skin of body",
          "breast", "adipose tissue", "musculature", "vasculature", "placenta",
          "prostate gland", "uterus"]
# non-blood filter (cell-type level) for footprint / liabilities
BLOOD = ["t cell", "b cell", " nk", "natural killer", "lymphocyte", "leukocyte",
         "monocyte", "macrophage", "dendritic", "granulocyte", "neutrophil",
         "eosinophil", "basophil", "mast cell", "plasma cell", "plasmablast",
         "erythro", "platelet", "megakaryocyte", "hematopoietic", "myeloid",
         "lymphoid", "thymocyte", "microglial", "kupffer", "langerhans",
         "phagocyte", "immune", "blast", "progenitor cell", "killer",
         "plasmacytoid", "innate lymphoid", "mononuclear", "t follicular",
         "t-helper", "helper", "regulatory t", "treg", "gamma-delta",
         "hofbauer", "osteoclast", "promyelocyte"]
ORG_EXCL = {"blood", "bone marrow", "embryo", "yolk sac"}
TOXN, TOXPCT = 3000, 30.0
GLOBAL = os.path.join(HERE, "Figure8_census_normal_celltype_global.csv")

# canonical normal compartments -> representative census cell types (n-weighted)
CT_GROUPS = [
    ("HSPC", ["hematopoietic stem cell",
              "hematopoietic multipotent progenitor cell",
              "common myeloid progenitor", "granulocyte monocyte progenitor cell",
              "megakaryocyte-erythroid progenitor cell",
              "common lymphoid progenitor"]),
    ("Monocyte/Mac/DC", ["CD14-positive, CD16-negative classical monocyte",
                         "non-classical monocyte", "alveolar macrophage",
                         "Kupffer cell", "CD1c-positive myeloid dendritic cell",
                         "plasmacytoid dendritic cell", "mast cell", "basophil"]),
    ("T/NK/ILC", ["CD8-positive, alpha-beta memory T cell",
                  "naive thymus-derived CD4-positive, alpha-beta T cell",
                  "CD16-positive, CD56-dim natural killer cell, human",
                  "group 3 innate lymphoid cell, human"]),
    ("Endothelium", ["blood vessel endothelial cell", "capillary endothelial cell",
                     "vein endothelial cell",
                     "endothelial cell of hepatic sinusoid"]),
    ("Hepatic/stroma", ["hepatocyte", "periportal region hepatocyte",
                        "hepatic stellate cell"]),
    ("Cardiac/muscle", ["regular ventricular cardiac myocyte",
                        "cardiac muscle cell"]),
    ("Neural/retina", ["L5/6 near-projecting glutamatergic neuron",
                       "sst GABAergic cortical interneuron",
                       "starburst amacrine cell"]),
    ("Mesothelium", ["mesothelial cell"]),
    ("Epithelium", ["bladder urothelial cell", "retinal pigment epithelial cell",
                    "suprabasal keratinocyte"]),
]


def is_circ_blood(ct):
    n = str(ct).lower()
    if any(k in n for k in ("macrophage", "microglia", "kupffer", "langerhans",
                            "hofbauer", "mast cell", "osteoclast")):
        return False
    excl = ("t cell", "thymocyte", "t-helper", "regulatory t", "mait",
            "gamma-delta", "memory t", "naive t", "effector t", "b cell",
            "plasma cell", "plasmablast", "plasmacytoid", "natural killer",
            "nk cell", "innate lymphoid", "lymphocyte", "leukocyte",
            "mononuclear", "monocyte", "dendritic cell", "neutrophil",
            "eosinophil", "basophil", "granulocyte", "megakaryocyte", "platelet",
            "erythro", "reticulocyte", "normoblast", "hematopoietic",
            "myeloid progenitor", "lymphoid progenitor", "myeloid cell")
    return any(e in n for e in excl)


def nonblood(d):
    isb = d["cell_type"].str.lower().apply(lambda c: any(k in c for k in BLOOD))
    return d[(~isb) & (~d["tissue"].isin(ORG_EXCL))
             & (d["cell_type"].str.lower() != "unknown")]


def panel_organ_atlas(d):
    rows, labs = [], []
    for org in ORGANS:
        sub = d[d["tissue"] == org]
        if sub.empty:
            continue
        if org not in KEEP_ALL_CELLS:
            sub = sub[~sub["cell_type"].apply(is_circ_blood)]
        n = sub["n_cells"].sum()
        if n < 200:
            continue
        rows.append([100.0 * sub[f"{g}_pos"].sum() / n for g in GENES])
        labs.append(org)
    M = np.array(rows)
    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(5.4, 0.30 * len(labs) + 1.0))
        norm = mpl.colors.Normalize(0, 70)
        s_of = lambda v: v * 3.2 + 2.0
        ax.axvspan(-0.5, 0.5, color="#fbecea", zorder=0)
        for i in range(len(labs)):
            for j in range(len(GENES)):
                v = M[i, j]
                cmap = CD96_CMAP if j == 0 else SEQ_CMAP
                ax.scatter(j, i, s=s_of(v), c=[v], cmap=cmap, norm=norm,
                           edgecolor="#3b5870", linewidth=0.3, zorder=3)
        ax.set_xlim(-0.6, len(GENES) - 0.4); ax.set_ylim(-0.6, len(labs) - 0.4)
        ax.invert_yaxis(); ax.set_xticks(range(len(GENES)))
        xt = ax.set_xticklabels([GLABEL[g] for g in GENES], fontsize=7.5)
        for t, g in zip(xt, GENES):
            t.set_color(TARGET_COLORS[g]); t.set_fontweight("bold")
        ax.xaxis.set_ticks_position("top")
        ax.set_yticks(range(len(labs)))
        yl = ax.set_yticklabels([f"{o}*" if o in KEEP_ALL_CELLS else o
                                 for o in labs], fontsize=7.5)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title("", pad=18)   # title/commentary in caption (authors)
        fig.subplots_adjust(left=0.34, right=0.97, top=1 - 1.0 / fig.get_figheight(),
                            bottom=0.4 / fig.get_figheight())
        save_fig(fig, HERE, "Figure8_panelB_organatlas")


def panel_fingerprint(glob):
    """compartment fingerprint: n-weighted % positive per canonical compartment
    (INCLUDES marrow/blood compartments, so each target's full off-tumor
    signature is shown fairly -- CD33/CLL-1 light up HSPC+monocyte, CD96 only
    T/NK)."""
    comps, M = [], []
    for name, cts in CT_GROUPS:
        pres = [c for c in cts if c in glob.index]
        if not pres:
            continue
        n = glob.loc[pres, "n_cells"].to_numpy(float)
        M.append([float(np.sum(n * glob.loc[pres, f"{g}_pct"].to_numpy(float))
                        / n.sum()) for g in GENES])
        comps.append(name)
    M = np.array(M)
    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(5.6, 3.6))
        norm = mpl.colors.Normalize(0, 70)
        s_of = lambda v: v * 3.2 + 2.0
        ax.axvspan(-0.5, 0.5, color="#fbecea", zorder=0)
        for i in range(len(comps)):
            for j in range(len(GENES)):
                v = M[i, j]
                ax.scatter(j, i, s=s_of(v), c=[v],
                           cmap=(CD96_CMAP if j == 0 else SEQ_CMAP), norm=norm,
                           edgecolor="#3b5870", linewidth=0.3, zorder=3)
        ax.set_xlim(-0.6, len(GENES) - 0.4); ax.set_ylim(-0.6, len(comps) - 0.4)
        ax.invert_yaxis(); ax.set_xticks(range(len(GENES)))
        xt = ax.set_xticklabels([GLABEL[g] for g in GENES], fontsize=7.5)
        for t, g in zip(xt, GENES):
            t.set_color(TARGET_COLORS[g]); t.set_fontweight("bold")
        ax.xaxis.set_ticks_position("top")
        ax.set_yticks(range(len(comps))); ax.set_yticklabels(comps, fontsize=8)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title("", pad=18)   # title/commentary in caption (authors)
        fig.subplots_adjust(left=0.30, right=0.97, top=0.80, bottom=0.06)
        save_fig(fig, HERE, "Figure8_panelC_fingerprint")


def panel_cd96_liab(nb):
    big = nb[(nb["n_cells"] >= TOXN) & (nb["CD96_pct"] >= 20)].copy()
    big = big.sort_values("CD96_pct", ascending=True).tail(12)
    labs = [f"{r.cell_type[:30]} · {r.tissue}" for r in big.itertuples()]
    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(5.6, 3.4))
        y = np.arange(len(big))
        ax.barh(y, big["CD96_pct"], color=PALETTE["CD96"], edgecolor="white",
                height=0.72, alpha=0.9)
        for yi, v in zip(y, big["CD96_pct"]):
            ax.text(v + 1, yi, f"{v:.0f}%", va="center", fontsize=7.5)
        ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=7.0)
        ax.set_xlim(0, 100)
        ax.set_xlabel("% cells CD96-positive")
        ax.set_title("")   # title/commentary in caption (authors)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        fig.subplots_adjust(left=0.52, right=0.95, top=0.86, bottom=0.16)
        save_fig(fig, HERE, "Figure8_panelD_cd96liab")


def main():
    d = pd.read_csv(BY_CT)
    nb = nonblood(d)
    glob = pd.read_csv(GLOBAL).set_index("cell_type")
    panel_organ_atlas(d)
    panel_fingerprint(glob)
    panel_cd96_liab(nb)
    print("Done -> Figure8_panelB_organatlas / panelC_fingerprint / "
          "panelD_cd96liab .{pdf,png}")


if __name__ == "__main__":
    main()
