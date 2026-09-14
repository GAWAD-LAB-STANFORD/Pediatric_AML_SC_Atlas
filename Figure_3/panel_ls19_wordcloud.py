#!/usr/bin/env python
"""Figure 3 (LS19) word cloud of ALL significantly enriched GO:BP pathways (GSEA, LS19 vs all other
leukemic cells; there are no significant down pathways, so every significant term is LS19-enriched).
Words come from the significant term names, sized by summed -log10(padj) across terms (dominant
themes largest) and coloured on the shared manuscript SEQ gradient (the dot-plot scheme): teal for
lesser weight -> warm/red for the strongest themes, so colour reinforces size. Data = LS19_gsea_sig.csv.
No fabrication."""
import os, sys, re, numpy as np, pandas as pd
from wordcloud import WordCloud
from matplotlib.colors import LinearSegmentedColormap, to_hex
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"); import config as C
OUT = os.path.join(C.BASE, "04_Main_Figures/Figure_3_panels")
sig = pd.read_csv(os.path.join(C.DATA, "LS19_gsea_sig.csv"))
STOP = {"of","to","the","in","a","and","or","by","via","process","positive","negative","regulation",
        "involved","cellular","response","establishment","maintenance","activity","mediated","pathway"}
weights = {}
for _, row in sig.iterrows():
    w = -np.log10(max(row["padj"], 1e-300))
    for tok in re.split(r"[^a-z0-9]+", str(row["term"]).lower()):
        if len(tok) > 2 and tok not in STOP:
            weights[tok] = weights.get(tok, 0.0) + w
# shared manuscript SEQ palette (dot-plot scheme); drop the two palest stops so words stay legible on white
SEQ = ["#001219","#005F73","#0A9396","#EE9B00","#CA6702","#AE2012","#9B2226"]
cmap = LinearSegmentedColormap.from_list("seq_wc", SEQ)
wmin, wmax = min(weights.values()), max(weights.values())
def seq_color(word, *a, **k):
    t = (weights.get(word, wmin) - wmin) / (wmax - wmin + 1e-9)
    return to_hex(cmap(0.15 + 0.85 * t))          # map weight -> SEQ position (bigger = warmer)
wc = WordCloud(width=1100, height=750, background_color="white", prefer_horizontal=0.9,
               color_func=seq_color, relative_scaling=0.5, max_words=60, random_state=1)
wc.generate_from_frequencies(weights)
fig, ax = plt.subplots(figsize=(6.2, 4.2)); ax.imshow(wc, interpolation="bilinear"); ax.axis("off")
for e in ("png", "pdf"): fig.savefig(os.path.join(OUT, f"Figure_3_LS19_wordcloud.{e}"), dpi=200, bbox_inches="tight")
print("wrote Figure_3_LS19_wordcloud from", len(sig), "significant UP terms; top words:",
      ", ".join(k for k, _ in sorted(weights.items(), key=lambda x: -x[1])[:12]))
