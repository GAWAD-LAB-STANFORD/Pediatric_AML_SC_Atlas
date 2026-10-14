# ====================================================================
# figS_ppac_headtohead.R  |  Figure S15 (per-PPAC efficacy vs toxicity)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : figS_ppac_headtohead.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript figS_ppac_headtohead.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Supplementary figure — per-PPAC head-to-head, EXACTLY like Figure 5C, six panels.
# Each panel = one PPAC group (All PPACs pooled + PPAC_1..5).  Same toxicity-vs-efficacy
# space as 5C, but efficacy is PPAC-specific:
#   x = leak     = % of normal HSPC / Myeloid_Pro positive (max)        [toxicity]
#   y = cov_pct  = % of that PPAC's cells positive                      [efficacy]
# size = rec_pct = % of contributing patients with >=20% coverage       [recurrence]
# Points = the TOP-5 PPAC-specific markers (by % of the PPAC's cells positive, from S17) for that
# PPAC + the 12 panel-E antigens (Figure 7B/C). Panel-E genes are coloured by group (leads identity,
# FLT3 control gold, clinical orange, PPAC-specific purple); discovered markers grey; CD96 crimson.
# Source: figS_ppac_headtohead.csv  (build_ppac_headtohead.py).  Like fig5.R panel C.
source("_style.R")

d <- sd("figS_ppac_headtohead.csv")
d$fill <- ifelse(d$label %in% names(TARGET12_COLORS),
                 as.character(TARGET12_COLORS[d$label]), PAL$grey)   # panel-E by group; discovered grey
GROUPS <- c("All PPAC", paste0("PPAC_", 1:5))
TITLES <- c("All PPACs (pooled)", paste0("PPAC_", 1:5))

panel <- function(g, title) {
  dd  <- d[d$group == g, ]
  flg <- isTRUE(dd$single_patient_flag[1])
  p <- ggplot(dd, aes(leak, cov_pct)) +
    annotate("rect", xmin = -2, xmax = 10, ymin = -Inf, ymax = Inf,
             fill = "#eafaf1", alpha = 0.7)            # sparing zone (passes toxicity funnel, leak<=10%)
  if (flg)
    p <- p + annotate("text", x = 56, y = 4, hjust = 1, vjust = 0, size = 2.4,
                      fontface = "italic", colour = "#a93226",
                      label = "single-patient PPAC\n(not generalizable)")
  p +
    geom_vline(xintercept = 10, linetype = "dotted", colour = "#7f8c8d", linewidth = 0.3) +
    geom_point(aes(size = rec_pct, fill = fill), shape = 21, colour = "black", stroke = 0.4) +
    scale_fill_identity() +
    scale_size_area(max_size = 7, limits = c(0, 100), name = "% patients\nrecurrent",
                    breaks = c(10, 25, 50, 75)) +
    ggrepel::geom_text_repel(aes(label = label, colour = fill,
                                 fontface = ifelse(label == "CD96", "bold", "plain")),
                             size = 2.5, max.overlaps = 30, seed = 1,
                             box.padding = 0.4, point.padding = 0.3,
                             min.segment.length = 0, segment.size = 0.2,
                             segment.colour = "grey65") +
    scale_colour_identity() +
    coord_cartesian(xlim = c(-2, 58), ylim = c(0, 85)) +
    labs(x = "% normal HSPC / Myeloid_Pro positive (max)  (toxicity)",
         y = "% of PPAC cells positive  (efficacy)", title = title,
         subtitle = sprintf("n = %s cells", format(dd$n_cells[1], big.mark = ","))) +
    theme_pub(8.5) +
    theme(plot.subtitle = element_text(size = 6.8, colour = "#7f8c8d"))
}

panels <- Map(panel, GROUPS, TITLES)

fig <- (panels[[1]] | panels[[2]]) /
       (panels[[3]] | panels[[4]]) /
       (panels[[5]] | panels[[6]]) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Per-PPAC efficacy vs toxicity - top-5 PPAC-specific markers + the 12 panel-E antigens",
    subtitle = paste0("Figure-5C construction applied within each poor-prognosis cluster. ",
                      "Each point is a surface marker: x = normal HSPC/Myeloid toxicity, ",
                      "y = % of that PPAC's cells covered, dot size = % of patients with recurrent ",
                      "(>=20%) coverage. Green = sparing zone (passes the toxicity funnel, leak <=10%). ",
                      "Panel-E antigens by group (clinical orange, etc.); discovered markers grey; CD96 crimson."),
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = ggtext::element_textbox_simple(size = 8, colour = "#34495e",
                                                                 margin = margin(b = 6)))) &
  theme(legend.position = "right")

savefig(fig, "supplementary/FigureS15_ppac_headtohead", 12, 13)
cat("Done -> figures/supplementary/FigureS15_ppac_headtohead.{pdf,png}\n")
