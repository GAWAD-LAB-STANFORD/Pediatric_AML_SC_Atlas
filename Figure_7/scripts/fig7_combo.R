# ====================================================================
# fig7_combo.R  |  Figure 7B (CD96 + one-partner additive window)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : fig7_combo.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript fig7_combo.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Figure 7 — additive coverage vs marrow toxicity of CD96 + ONE partner ("old 6C" format).
# Two Pareto panels: WITHIN patients (% blasts covered, mean over the 19 CD96-targetable pts)
# and BETWEEN patients (% pts with >=20% blasts covered). x = % normal HSPC/Myeloid-Pro hit by
# (CD96 OR partner). Green zone = marrow-sparing (<10%), red zone = marrow-toxic. CD96 = base
# (crimson diamond); each partner is one addition. Red ring = highly toxic on a vital organ in
# Figure 6 (CD9, CD123, IL1RAP, ABCA7) -> disqualified despite where they land here.
source("_style.R")
df <- as.data.frame(sd("fig7_combo.csv"))
df$group[df$label == "CD7"] <- "clinical"          # corrected scheme: CD7 = clinical CAR-T antigen
base <- df[df$base, ]; part <- df[!df$base, ]
bx <- base$tox[1]; TOXGATE <- 10

GCOL <- c(lead = "#2980b9", clinical = "#e8590c", control = "#b7950b", ppac = "#7f8c8d")
GLAB <- c(lead = "lead candidate", clinical = "clinical antigen",
          control = "FLT3 (pan-AML control)", ppac = "additional candidate")
part$group <- factor(part$group, levels = names(GCOL))

mkpanel <- function(ycol, ylab, ttl) {
  by <- base[[ycol]][1]
  yr <- range(df[[ycol]]); pad <- diff(yr) * 0.16
  ylo <- floor((yr[1] - pad) / 5) * 5; yhi <- ceiling((yr[2] + pad) / 5) * 5
  ph <- part[part$highly_toxic, ]
  seg <- data.frame(x = bx, y = by, xend = part$tox, yend = part[[ycol]])
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = TOXGATE, ymin = -Inf, ymax = Inf, fill = "#eafaf1") +
    annotate("rect", xmin = TOXGATE, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#fdecea") +
    geom_vline(xintercept = TOXGATE, linetype = "dashed", colour = "#c0392b", linewidth = 0.4) +
    annotate("text", x = TOXGATE - 0.6, y = yhi, label = "marrow-sparing", hjust = 1, vjust = 1,
             size = 2.7, colour = "#1e8449", fontface = "italic") +
    annotate("text", x = TOXGATE + 0.6, y = yhi, label = "marrow toxicity", hjust = 0, vjust = 1,
             size = 2.7, colour = "#c0392b", fontface = "italic") +
    # each partner is an addition to CD96 -> faint fan from the base point
    geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "#c2c9cf", linewidth = 0.25, linetype = "22") +
    # red ring = highly toxic on a vital organ (Figure 6)
    geom_point(data = ph, aes(tox, .data[[ycol]]), shape = 21, fill = NA,
               colour = "#c0392b", size = 7.2, stroke = 1.35) +
    geom_point(data = part, aes(tox, .data[[ycol]], fill = group), shape = 21,
               size = 3.7, stroke = 0.4, colour = "white") +
    geom_point(data = base, aes(tox, .data[[ycol]]), shape = 23, fill = "#c0392b",
               colour = "white", size = 5, stroke = 0.5) +
    ggrepel::geom_text_repel(data = part, aes(tox, .data[[ycol]], label = label, colour = group),
                             size = 3, fontface = "bold", seed = 7, box.padding = 0.55,
                             point.padding = 0.4, max.overlaps = Inf, min.segment.length = 0,
                             segment.size = 0.2, segment.colour = "#aab2b8", show.legend = FALSE) +
    ggrepel::geom_text_repel(data = base, aes(tox, .data[[ycol]], label = label), colour = "#c0392b",
                             size = 3.15, fontface = "bold", nudge_y = -(yhi - ylo) * 0.05,
                             nudge_x = 1.6, seed = 7, segment.colour = NA) +
    scale_fill_manual(values = GCOL, labels = GLAB, name = NULL, drop = FALSE) +
    scale_colour_manual(values = GCOL, guide = "none") +
    scale_x_continuous(breaks = seq(0, 40, 10)) +
    coord_cartesian(xlim = c(0, 39), ylim = c(ylo, yhi), clip = "off") +
    labs(title = ttl, y = ylab,
         x = "% normal HSPC / Myeloid-Pro positive  ·  CD96 OR partner  (marrow toxicity →)") +
    theme_pub(10) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 11.5))
}

pL <- mkpanel("within", "% leukemic blasts covered by CD96 + partner\n(mean across 19 CD96-targetable patients · single cell)",
              "Within patients")
pR <- mkpanel("between", "% patients targetable by CD96 + partner\n(≥20% of blasts covered · single cell)",
              "Between patients")

fig <- (pL | pR) +
  plot_layout(guides = "collect") &           # no title/subtitle; commentary -> caption
  theme(legend.position = "bottom")

savefig(fig, "Figure7_CD96_additive_window", 13.5, 7)
cat("\nCD96 + partner (sorted by within-patient depth):\n")
print(df[order(-df$within), c("label", "group", "within", "between", "tox", "highly_toxic")], row.names = FALSE)
