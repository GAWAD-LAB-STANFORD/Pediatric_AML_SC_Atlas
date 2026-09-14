# ====================================================================
# figS_ls_headtohead.R  |  Per-poor-prognosis-STATE (LS) efficacy vs toxicity
# Replaces the retired per-PPAC version (figS_ppac_headtohead.R) on the 49-state scheme.
# Figure-5C construction applied within each poor-prognosis LEUKEMIC STATE (LS), six panels:
#   x = leak    = % normal HSPC/Myeloid_Pro positive (max)   [toxicity]
#   y = cov_pct = % of that state's cells positive           [efficacy]
# size = rec_pct = % of contributing patients with >=20% coverage [recurrence]
# Points = the top-5 state-specific markers (from funnel_ls_specific.csv) + the 12 panel-E
# antigens. Panel-E genes coloured by group (TARGET12_COLORS); discovered grey; CD96 crimson.
# Source: figS_ls_headtohead.csv  (build_ls_headtohead.py). Mirrors fig5.R panel C.
# ====================================================================
source("_style.R")

d <- sd("figS_ls_headtohead.csv")
d$fill <- ifelse(d$label %in% names(TARGET12_COLORS),
                 as.character(TARGET12_COLORS[d$label]), PAL$grey)
GROUPS <- c("All poor LS", "LS_10", "LS_19", "LS_20", "LS_21", "LS_27")
TITLES <- c("All poor-prognosis LS (pooled)", "LS_10 (proliferative)", "LS_19 (pDC)",
            "LS_20 (proliferative)", "LS_21 (macrophage)", "LS_27 (GMP)")

panel <- function(g, title) {
  dd  <- d[d$group == g, ]
  flg <- isTRUE(dd$single_patient_flag[1])
  p <- ggplot(dd, aes(leak, cov_pct)) +
    annotate("rect", xmin = -2, xmax = 10, ymin = -Inf, ymax = Inf, fill = "#eafaf1", alpha = 0.7)
  if (flg)
    p <- p + annotate("text", x = 56, y = 4, hjust = 1, vjust = 0, size = 2.4,
                      fontface = "italic", colour = "#a93226",
                      label = "single-patient state\n(not generalizable)")
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
                             min.segment.length = 0, segment.size = 0.2, segment.colour = "grey65") +
    scale_colour_identity() +
    coord_cartesian(xlim = c(-2, 58), ylim = c(0, 85)) +
    labs(x = "% normal HSPC / Myeloid_Pro positive (max)  (toxicity)",
         y = "% of state cells positive  (efficacy)", title = title,
         subtitle = sprintf("n = %s cells", format(dd$n_cells[1], big.mark = ","))) +
    theme_pub(8.5) +
    theme(plot.subtitle = element_text(size = 6.8, colour = "#7f8c8d"))
}

panels <- Map(panel, GROUPS, TITLES)
fig <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) / (panels[[5]] | panels[[6]]) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Per-state efficacy vs toxicity - top-5 state-specific markers + the 12 panel-E antigens",
    subtitle = NULL,
    theme = theme(plot.title = element_text(face = "bold", size = 12),
                  plot.subtitle = ggtext::element_textbox_simple(size = 8, colour = "#34495e",
                                                                 margin = margin(b = 6)))) &
  theme(legend.position = "right")

savefig(fig, "supplementary/FigureS16_ls_headtohead", 12, 13)
cat("Done -> figures/supplementary/FigureS16_ls_headtohead.{pdf,png}\n")
