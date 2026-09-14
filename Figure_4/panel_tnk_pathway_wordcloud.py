#!/usr/bin/env python
"""Figure 4: word cloud of the GO:BP pathways enriched in AML vs healthy-BM T/NK cells (GSEA on the full
gene ranking; 229 significant UP pathways, so a word cloud summarises the dominant themes better than a
bar chart). Words from the significant UP term names, sized by summed -log10(padj), coloured on the shared
SEQ gradient. Data = tnk_gsea_all.csv. No fabrication."""
import os, re, numpy as np, pandas as pd
from wordcloud import WordCloud
from matplotlib.colors import LinearSegmentedColormap, to_hex
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
V2  = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT = "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
g = pd.read_csv(os.path.join(V2, "tnk_gsea_all.csv"))
sig = g[(g["padj"] < 0.05) & (g["NES"] > 0)].copy()
STOP = {"of","to","the","in","a","and","or","by","via","process","positive","negative","regulation",
        "involved","cellular","response","establishment","maintenance","activity","mediated","pathway",
        "organization"}
weights = {}
for _, row in sig.iterrows():
    w = -np.log10(max(row["padj"], 1e-300))
    for tok in re.split(r"[^a-z0-9]+", str(row["term"]).lower()):
        if len(tok) > 2 and tok not in STOP:
            weights[tok] = weights.get(tok, 0.0) + w
SEQ = ["#001219","#005F73","#0A9396","#EE9B00","#CA6702","#AE2012","#9B2226"]
cmap = LinearSegmentedColormap.from_list("seq_wc", SEQ)
wmin, wmax = min(weights.values()), max(weights.values())
def seq_color(word, *a, **k):
    t = (weights.get(word, wmin) - wmin) / (wmax - wmin + 1e-9)
    return to_hex(cmap(0.15 + 0.85 * t))
wc = WordCloud(width=1100, height=750, background_color="white", prefer_horizontal=0.9,
               color_func=seq_color, relative_scaling=0.5, max_words=70, random_state=1)
wc.generate_from_frequencies(weights)
fig, ax = plt.subplots(figsize=(6.4, 4.4)); ax.imshow(wc, interpolation="bilinear"); ax.axis("off")
for e in ("png", "pdf"): fig.savefig(os.path.join(OUT, f"Figure_4_tnk_pathway_wordcloud.{e}"), dpi=200, bbox_inches="tight")
print("wrote Figure_4_tnk_pathway_wordcloud from", len(sig), "significant UP pathways; top words:",
      ", ".join(k for k, _ in sorted(weights.items(), key=lambda x: -x[1])[:14]))
