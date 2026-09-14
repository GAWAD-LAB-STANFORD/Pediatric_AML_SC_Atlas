# ====================================================================
# fig7_allpairs.R  |  Figure 7D-E & Figure S26 (all 66 antigen pairs)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : fig7_allpairs.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript fig7_allpairs.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Every 2-target combination of the 12 panel-E antigens (66 pairs) — toxicity (x) vs efficacy (y).
# NOT anchored on CD96. Winners = Pareto frontier (best efficacy reachable at each toxicity),
# filled by marrow toxicity and labelled with BOTH partners.  Renders, for the deck:
#   Figure7D_punchline = consolidated landscape (scatter); Figure7E_winners = winners QUANTIFIED
#   (ranked bars of the balanced-frontier combos).  And the supplement S26 (four views below).
# Three scatter views:
#   A  within patients (depth)   = % blasts covered, mean across 28 patients
#   B  between patients (breadth) = % patients with >=20% of blasts covered
#   C  CONSOLIDATED balanced score = 100 x sqrt( depth_norm x breadth_norm ), each normalised to
#      its max over pairs — high only when a combo is strong on BOTH axes (=> the best overall combo).
# Low-toxicity (<10%) winners are detailed with (depth / breadth / toxicity %). Green zone = sparing;
# dark-red ring = a winner with a TOXIC member: organ-toxic in Figure 6 (CD9/CD123/IL1RAP/ABCA7)
# OR the marrow-toxic pan-AML control FLT3 (so every FLT3 combination is ringed).
source("_style.R")
d <- as.data.frame(sd("fig7_allpairs.csv"))
TOXIC <- c("CD9", "CD123", "IL1RAP", "ABCA7", "FLT3")    # organ-toxic (Fig 6) + FLT3 marrow-toxic control
d$tox_flag <- d$g1 %in% TOXIC | d$g2 %in% TOXIC
pairs <- d[d$kind == "pair", ]
mw <- max(pairs$within); mb <- max(pairs$across)
d$balanced <- round(100 * sqrt((d$within / mw) * (d$across / mb)), 1)
pairs <- d[d$kind == "pair", ]; singles <- d[d$kind == "single", ]
TOXGATE <- 10; TOXMAX <- ceiling(max(pairs$toxicity))

# one shared marrow-toxicity colour scale for every panel of the all-pairs figure
tox_fill <- function() scale_fill_gradientn(
  colours = c("#1e8449", "#7dc47f", "#f7f7bf", "#ef9b5a", "#c0392b"),
  values = scales::rescale(c(0, 6, 10, 25, TOXMAX)), limits = c(0, TOXMAX),
  name = "marrow toxicity\n(% HSPC/Myeloid⁺)")

mkpanel <- function(ycol, ylab, ttl) {
  fr <- pairs[order(pairs$toxicity), ]; fr <- fr[fr[[ycol]] >= cummax(fr[[ycol]]) - 1e-9, ]   # winners
  other <- pairs[!(pairs$combo %in% fr$combo), ]
  fr$lab <- ifelse(fr$toxicity < TOXGATE,                       # low-tox winners get the numbers
                   sprintf("%s\n(%.0f / %.0f / %.0f)", fr$combo, fr$within, fr$across, fr$toxicity),
                   fr$combo)
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = TOXGATE, ymin = -Inf, ymax = Inf, fill = "#eafaf1") +
    annotate("rect", xmin = TOXGATE, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#fdecea") +
    geom_vline(xintercept = TOXGATE, linetype = "dashed", colour = "#c0392b", linewidth = 0.35) +
    geom_line(data = fr[order(fr$toxicity), ], aes(toxicity, .data[[ycol]]),
              colour = "#7f8c8d", linetype = "dotted", linewidth = 0.5) +
    geom_point(data = singles, aes(toxicity, .data[[ycol]]), shape = 23, fill = NA,
               colour = "#b0b8bf", size = 1.9, stroke = 0.35) +
    geom_point(data = other, aes(toxicity, .data[[ycol]]), shape = 21, fill = "#d6dbe1",
               colour = "grey62", size = 2.4, stroke = 0.2) +
    geom_point(data = fr[fr$tox_flag, ], aes(toxicity, .data[[ycol]]), shape = 21, fill = NA,
               colour = "#7b241c", size = 7.0, stroke = 1.3) +
    geom_point(data = fr, aes(toxicity, .data[[ycol]], fill = toxicity), shape = 21,
               colour = "grey20", size = 4.3, stroke = 0.4) +
    tox_fill() +
    ggrepel::geom_text_repel(data = fr, aes(toxicity, .data[[ycol]], label = lab),
                             size = 2.45, fontface = "bold", colour = "#1c2833", lineheight = 0.86,
                             max.overlaps = Inf, box.padding = 0.5, point.padding = 0.4, seed = 4,
                             min.segment.length = 0, segment.size = 0.2, segment.colour = "#aab2b8") +
    coord_cartesian(ylim = c(0, NA)) +
    labs(x = "normal marrow toxicity  (% HSPC / Myeloid-Pro positive, X OR Y)", y = ylab, title = ttl) +
    theme_pub(10) + theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 11))
}
# ---- winners quantified, SEPARATED by toxicity tier ------------------------------------------
# Two faceted tiers ranked by balanced coverage score:
#   top    = LOW-TOXICITY options (organ-clean AND marrow-sparing <10%)  -> the deployable combos
#   bottom = HIGHER-COVERAGE but TOXIC (organ-toxic/FLT3 member OR marrow-suppressive >=10%) -> ceiling
# Bar fill = marrow toxicity; tip = score · depth / breadth / tox%, printed dark-red when the pair
# carries an organ-toxic (Fig 6) or FLT3 member, dark-green when clean.
LTOP <- "Low-toxicity options    ·    organ-clean & marrow-sparing (<10%)"
HTOP <- "Higher coverage    ·    toxic member (organ-toxic / FLT3)  or  marrow-suppressive (≥10%)"
mkbar <- function(ttl, ntier = 6) {
  cl  <- pairs[!pairs$tox_flag & pairs$toxicity < TOXGATE, ]   # organ-clean AND marrow-sparing
  lo  <- head(cl[order(-cl$balanced), ], ntier);  lo$tier <- LTOP
  fr  <- pairs[order(pairs$toxicity), ]
  fr  <- fr[fr$balanced >= cummax(fr$balanced) - 1e-9, ]       # Pareto frontier
  hi  <- fr[fr$tox_flag | fr$toxicity >= TOXGATE, ]
  hi  <- head(hi[order(-hi$balanced), ], ntier); hi$tier <- HTOP
  w   <- rbind(lo, hi); w$tier <- factor(w$tier, levels = c(LTOP, HTOP))
  w   <- w[order(w$tier, w$balanced), ]; w$combo <- factor(w$combo, levels = w$combo)
  w$tip <- sprintf("%.0f  ·  %.0f / %.0f / %.0f%%", w$balanced, w$within, w$across, w$toxicity)
  ggplot(w, aes(balanced, combo)) +
    geom_col(aes(fill = toxicity), width = 0.7, colour = "grey25", linewidth = 0.3) +
    geom_text(aes(label = tip, colour = tox_flag), hjust = -0.06, size = 2.7, fontface = "bold") +
    tox_fill() +
    scale_colour_manual(values = c(`TRUE` = "#7b241c", `FALSE` = "#145a32"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.38)), limits = c(0, NA)) +
    facet_wrap(~ tier, ncol = 1, scales = "free_y") +
    labs(x = "balanced coverage score   (bar fill = marrow toxicity %;  tip = score · depth / breadth / toxicity %)",
         y = NULL, title = ttl) +
    theme_pub(10) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
          axis.text.y = element_text(face = "bold", size = 9),
          strip.text = element_text(face = "bold", size = 10, hjust = 0),
          strip.background = element_rect(fill = "grey92", colour = NA),
          plot.title = element_text(face = "bold", size = 11))
}

pW <- mkpanel("within",   "% leukemic blasts covered (X OR Y)\n(mean across 28 AML patients)", "A   Within patients (depth)")
pB <- mkpanel("across",   "% patients targetable (X OR Y)\n(≥20% of blasts covered)", "B   Between patients (breadth)")
pC <- mkpanel("balanced", "balanced coverage score\n(geom. mean of depth × breadth, 0–100)",
              "C   Consolidated — best combo on both axes   (low-tox labels show depth / breadth / toxicity %)")
pD <- mkbar("D   Winners quantified — low-toxicity options vs higher-coverage toxic ceiling")

fig <- (pW | pB) / (pC | pD) + plot_layout(heights = c(1, 1.18), guides = "collect") &
  theme(legend.position = "bottom")
savefig(fig, "supplementary/FigureS26_allpairs_window", 16, 12.8)

# ---- Figure 7 panel D (punch line): consolidated landscape, title-less for the assembler -------
pDsolo <- mkpanel("balanced", "balanced coverage score\n(geom. mean of depth × breadth, 0–100)", NULL) +
  theme(legend.position = "bottom")
savefig(pDsolo, "Figure7D_punchline", 15.5, 6)

# ---- Figure 7 panel E (winners quantified, tier-separated): ranked bars, title-less ------------
pEsolo <- mkbar(NULL) + theme(legend.position = "bottom")
savefig(pEsolo, "Figure7E_winners", 16, 6.0)

cat("\nConsolidated balanced-score ranking (top 12; depth/breadth/tox):\n")
print(head(pairs[order(-pairs$balanced), c("combo", "within", "across", "toxicity", "balanced", "tox_flag")], 12),
      row.names = FALSE)
cat("\nBest MARROW-SPARING (<10%) combos by balanced score:\n")
ms <- pairs[pairs$toxicity < TOXGATE, ]
print(head(ms[order(-ms$balanced), c("combo", "within", "across", "toxicity", "balanced", "tox_flag")], 8),
      row.names = FALSE)
