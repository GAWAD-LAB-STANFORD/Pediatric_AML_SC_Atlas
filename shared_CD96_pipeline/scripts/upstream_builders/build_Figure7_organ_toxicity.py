#!/usr/bin/env python
# ====================================================================
# build_Figure7_organ_toxicity.py  |  Upstream builder (regenerates source_data from raw inputs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : python build_Figure7_organ_toxicity.py   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

"""Body-wide OFF-TUMOR toxicity map across NON-blood tissues, separated by organ.

Blood/marrow are covered elsewhere (Fig 6, the differentiation panel), so here we
drop all blood/immune cell types and blood/marrow organs and focus on the major
parenchymal/structural cell types throughout the body where AT LEAST ONE of the 9
targets is predicted toxic (>= TOX% of cells positive).  Rows are grouped and
separated BY ORGAN; columns are the 9 targets (CD96 in its red colormap).

Filters: non-blood cell type; organ not blood/marrow/embryo; n_cells >= MINN
(major population); max marker >= TOX; top TOPN cell types per organ.

Message: the body-wide toxicity is driven by GPR56 / CD123 / IL1RAP / MSLN /
FLT3; CD96 is clean across nearly every organ (only urothelium + adipocytes).
Source: CZ CELLxGENE Census per-(tissue,cell-type) table.  No fabricated values.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.transforms import blended_transform_factory as blend

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, "supplements")
sys.path.insert(0, SUPP)
from _fig7_style import PALETTE, RC, save_fig, SEQ_CMAP, CD96_CMAP, TARGET_COLORS

BY_CT = os.path.join(HERE, "Figure8_census_normal_by_celltype.csv")
GLOBAL_CT = os.path.join(HERE, "Figure8_census_normal_celltype_global.csv")
HEME = [("hematopoietic stem cell", "HSC"),
        ("hematopoietic multipotent progenitor cell", "MPP"),
        ("common myeloid progenitor", "CMP"),
        ("granulocyte monocyte progenitor cell", "GMP"),
        ("megakaryocyte-erythroid progenitor cell", "MEP"),
        ("common lymphoid progenitor", "CLP")]
GENES = ["CD96", "CD33", "IL3RA", "CLEC12A", "FLT3", "MSLN", "CD70", "ADGRG1",
         "IL1RAP"]
GLABEL = {"CD96": "CD96", "CD33": "CD33", "IL3RA": "CD123", "CLEC12A": "CLL-1",
          "FLT3": "FLT3", "MSLN": "MSLN", "CD70": "CD70", "ADGRG1": "GPR56",
          "IL1RAP": "IL1RAP"}
MINN, TOX, TOPN, VMAX = 3000, 30.0, 2, 80.0

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

ORGAN_SHORT = {"Bone marrow (HSC→prog.)": "BONE MARROW\n(HSC→prog.)",
               "central nervous system": "CNS", "respiratory system": "resp.",
               "reproductive system": "reprod.", "digestive system": "digest.",
               "urinary bladder": "bladder", "skin of body": "skin",
               "large intestine": "lg intestine", "small intestine": "sm intestine"}


def main():
    d = pd.read_csv(BY_CT)
    d["maxpct"] = d[[f"{g}_pct" for g in GENES]].max(axis=1)
    is_blood = d["cell_type"].str.lower().apply(
        lambda c: any(k in c for k in BLOOD))
    nb = d[(~is_blood) & (~d["tissue"].isin(ORG_EXCL))
           & (d["cell_type"].str.lower() != "unknown")
           & (d["n_cells"] >= MINN) & (d["maxpct"] >= TOX)].copy()
    keep = (nb.sort_values("maxpct", ascending=False)
              .groupby("tissue").head(TOPN))
    # order: organ alphabetical, within organ by max marker desc
    keep = keep.sort_values(["tissue", "maxpct"], ascending=[True, False])

    # prepend the haematopoietic differentiation stages (HSC -> progenitors,
    # incl. HSC & MPP) as a top "bone marrow" group, in biological order, so the
    # stem/progenitor toxicity stays in view alongside the body-wide tissues
    glob = pd.read_csv(GLOBAL_CT).set_index("cell_type")
    hrows = []
    for ct, lab in HEME:
        r = {"tissue": "Bone marrow (HSC→prog.)", "cell_type": lab,
             "n_cells": int(glob.loc[ct, "n_cells"])}
        for g in GENES:
            r[f"{g}_pct"] = float(glob.loc[ct, f"{g}_pct"])
            r[f"{g}_mean"] = float(glob.loc[ct, f"{g}_mean"])
        r["maxpct"] = max(r[f"{g}_pct"] for g in GENES)
        hrows.append(r)
    keep = pd.concat([pd.DataFrame(hrows), keep], ignore_index=True, sort=False)
    keep.to_csv(os.path.join(HERE, "Figure7_organ_toxicity.csv"), index=False)
    nrow = len(keep)
    print(f"  {nrow} cell types across {keep['tissue'].nunique()} organs "
          f"(non-blood, n≥{int(MINN)}, max≥{int(TOX)}%, top{TOPN}/organ)")

    M = keep[[f"{g}_pct" for g in GENES]].to_numpy()
    cts = keep["cell_type"].tolist()
    orgs = keep["tissue"].tolist()
    ns = keep["n_cells"].tolist()

    cmap = SEQ_CMAP
    norm = mpl.colors.Normalize(vmin=0, vmax=VMAX)

    def s_of(v):
        return v * 3.4 + 2.0

    with mpl.rc_context(RC):
        fig, ax = plt.subplots(figsize=(9.8, 0.205 * nrow + 1.6))
        ax.set_axisbelow(True)
        ax.grid(color="#eef2f6", lw=0.5, zorder=1)
        # CD96 column band
        ax.axvspan(-0.5, 0.5, color="#fbecea", zorder=0)
        xs, ys, ss, cc = [], [], [], []
        rx, ry, rs, rc = [], [], [], []
        for i in range(nrow):
            for j, gn in enumerate(GENES):
                v = M[i, j]
                if j == 0:                      # CD96 -> red
                    rx.append(j); ry.append(i); rs.append(s_of(v)); rc.append(v)
                else:
                    xs.append(j); ys.append(i); ss.append(s_of(v)); cc.append(v)
        ax.scatter(xs, ys, s=ss, c=cc, cmap=cmap, norm=norm, edgecolor="#3b5870",
                   linewidth=0.3, zorder=3)
        ax.scatter(rx, ry, s=rs, c=rc, cmap=CD96_CMAP, norm=norm,
                   edgecolor="#7d2018", linewidth=0.4, zorder=4)

        ax.set_xlim(-0.6, len(GENES) - 0.4); ax.set_ylim(-0.6, nrow - 0.4)
        ax.invert_yaxis()
        ax.set_xticks(range(len(GENES)))
        xt = ax.set_xticklabels([GLABEL[g] for g in GENES], fontsize=8.5)
        for t, g in zip(xt, GENES):
            t.set_color(TARGET_COLORS.get(g, "#333")); t.set_fontweight("bold")
        ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
        ax.set_yticks(range(nrow))
        ax.set_yticklabels([f"{c[:34]}  ({n/1000:.0f}k)" for c, n in zip(cts, ns)],
                           fontsize=6.3)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)

        # organ group separators + left-margin organ labels
        tf = blend(ax.transAxes, ax.transData)
        start = 0
        for k in range(1, nrow + 1):
            if k == nrow or orgs[k] != orgs[start]:
                if start != 0:
                    ax.axhline(start - 0.5, color="#c9d3dd", lw=0.7, zorder=2)
                mid = (start + k - 1) / 2.0
                name = ORGAN_SHORT.get(orgs[start], orgs[start])
                ax.text(-0.46, mid, name, transform=tf, ha="left", va="center",
                        fontsize=7.0, fontweight="bold", color="#1c2833",
                        clip_on=False)
                start = k

        # size+colour legend
        axleg = fig.add_axes([0.78, 0.005, 0.20, 0.04]); axleg.axis("off")
        for v in (20, 50, 80):
            axleg.scatter([], [], s=s_of(v), color=cmap(norm(v)),
                          edgecolor="#3b5870", linewidth=0.4, label=f"{v}%")
        axleg.legend(title="% cells positive", loc="center", ncol=3,
                     frameon=False, fontsize=7.5, title_fontsize=7.5,
                     columnspacing=1.4, handletextpad=0.6)

        ax.set_title("", pad=18)   # title/commentary in caption (authors)
        fig.subplots_adjust(left=0.38, right=0.985, top=1 - 0.6 / fig.get_figheight(),
                            bottom=0.6 / fig.get_figheight())
        save_fig(fig, HERE, "Figure7_organ_toxicity")
    print("Done -> Figure7_organ_toxicity.{pdf,png,csv}")


if __name__ == "__main__":
    main()
