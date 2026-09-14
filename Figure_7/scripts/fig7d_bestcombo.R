# ====================================================================
# fig7d_bestcombo.R  |  Figure 7C (best second target decision view)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : fig7_combo.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript fig7d_bestcombo.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Figure 7D — choosing the best second target for CD96 by integrating WITHIN-patient depth and
# ACROSS-patient breadth into one decision view. Each point = CD96 + one partner.
#   x = across-patient breadth (% patients targetable, >=20% blasts covered by CD96 OR partner)
#   y = within-patient depth   (% blasts covered, mean across the 19 CD96-targetable patients)
#   fill  = marrow toxicity (% normal HSPC/Myeloid_Pro hit) ; red ring = organ-toxic in Figure 6
#   green dashes = best toxicity-free Pareto frontier (marrow-sparing AND organ-clean partners)
# Best combo = top-right (high depth + breadth), green fill, no ring.
source("_style.R")
df <- as.data.frame(sd("fig7_combo.csv"))
df$gene <- sub("^\\+", "", df$label)
ORGAN_TOXIC <- c("CD9", "CD123", "IL1RAP", "ABCA7")     # organ-toxic disqualified set, Figure 6
df$organ <- df$gene %in% ORGAN_TOXIC
base <- df[df$base, ]; part <- df[!df$base, ]
bx <- base$between[1]; by <- base$within[1]

# Two best 2-target combinations per strategy — each line = CD96 + ONE partner (no multi-target paths):
#   green = limit toxicity  (CD96+SUCNR1, CD96+ITGAX)
#   red   = maximize efficacy (CD96+CD33, CD96+CLL-1)
gseg <- transform(part[part$gene %in% c("SUCNR1", "ITGAX"), ], x0 = bx, y0 = by)
rseg <- transform(part[part$gene %in% c("CD33", "CLL-1"), ], x0 = bx, y0 = by)

p <- ggplot() +
  geom_segment(data = rseg, aes(x = x0, y = y0, xend = between, yend = within),
               colour = "#c0392b", linewidth = 0.7, linetype = "22") +
  geom_segment(data = gseg, aes(x = x0, y = y0, xend = between, yend = within),
               colour = "#1e8449", linewidth = 0.7, linetype = "22") +
  geom_point(data = part[part$organ, ], aes(between, within), shape = 21, fill = NA,
             colour = "#c0392b", size = 8.7, stroke = 1.4) +
  geom_point(data = part, aes(between, within, fill = tox), shape = 21,
             size = 5.3, stroke = 0.4, colour = "white") +
  geom_point(data = base, aes(between, within), shape = 23, fill = "#c0392b",
             size = 5.6, stroke = 0.5, colour = "white") +
  scale_fill_gradientn(colours = c("#1e8449", "#7dc47f", "#f7f7bf", "#ef9b5a", "#c0392b"),
                       values = scales::rescale(c(0, 6, 10, 20, 37)), limits = c(0, 37),
                       name = "marrow toxicity\n(% HSPC/Myeloid⁺)") +
  ggrepel::geom_text_repel(data = df, aes(between, within, label = label,
                           fontface = ifelse(base, "bold.italic", "bold")), size = 3.1,
                           box.padding = 0.55, point.padding = 0.4, max.overlaps = Inf, seed = 3,
                           min.segment.length = 0, segment.colour = "#aab2b8", segment.size = 0.2) +
  annotate("text", x = 83.5, y = 42.8, label = "① limit toxicity\n+ITGAX (breadth)\n+SUCNR1 (depth)",
           colour = "#1e8449", size = 2.95, fontface = "bold", hjust = 0, lineheight = 0.9) +
  annotate("text", x = 76.5, y = 64.2, label = "② maximize efficacy\n+CD33 / +CLL-1\n(accepts marrow toxicity)",
           colour = "#c0392b", size = 2.95, fontface = "bold", hjust = 0, lineheight = 0.9) +
  scale_x_continuous(limits = c(58, 97), breaks = seq(60, 95, 10)) +
  scale_y_continuous(limits = c(38, 70), breaks = seq(40, 70, 10)) +
  labs(x = "Across patients →  % targetable by CD96 + partner  (≥20 % blasts)",
       y = "Within patients →  % blasts covered by CD96 + partner\n(mean across 19 CD96-targetable patients)") +
  theme_pub(11) +                              # no title/subtitle; commentary -> caption
  theme(legend.position = c(0.012, 0.985), legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", colour = "grey85", linewidth = 0.3),
        legend.key.height = unit(0.8, "cm"), legend.key.width = unit(0.3, "cm"),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7))

savefig(p, "Figure7D_best_combo", 9.5, 7.8)
cat("\nCD96 + partner, ranked by within+across coverage (organ-toxic flagged):\n")
print(df[order(-(df$within + df$between)), c("gene", "group", "within", "between", "tox", "organ")], row.names = FALSE)
