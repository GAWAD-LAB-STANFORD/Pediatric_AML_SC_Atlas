# ====================================================================
# figS_ls_window.R  |  Additive therapeutic window per poor-prognosis STATE (LS)
# Replaces AND MERGES the retired per-PAC (figS_pac_window.R) + per-PPAC (figS_ppac_additive.R)
# versions on the corrected 49-state scheme (S27 dropped as a duplicate). Sparing zone <=10% (Fig 7).
# Each point = CD96 alone or CD96 OR one partner; x = normal HSPC/Myeloid toxicity,
# y = % of that state's cells covered (efficacy); dotted = Pareto frontier (go by toxicity,
# read the efficacy reached). A = all poor LS pooled; per-state facets = the 5 multi-patient
# poor states. Source: figS_ls_window.csv (build_ls_window.py). No figure subtitle (see caption).
# ====================================================================
source("_style.R")

d <- sd("figS_ls_window.csv")
d$label <- factor(d$label, levels = names(TARGET12_COLORS))
GREEN <- "#eafaf1"; RED <- "#fdecea"
zone <- list(annotate("rect", xmin = -4, xmax = 10, ymin = -Inf, ymax = Inf, fill = GREEN),
             annotate("rect", xmin = 10, xmax = 115, ymin = -Inf, ymax = Inf, fill = RED))
frontier <- function(dd) { dd <- dd[order(dd$toxicity), ]; dd[dd$efficacy >= cummax(dd$efficacy) - 1e-9, ] }
PROG <- c(LS_10 = "LS_10 (proliferative)", LS_19 = "LS_19 (pDC)", LS_20 = "LS_20 (proliferative)",
          LS_21 = "LS_21 (macrophage)", LS_27 = "LS_27 (GMP)")

pA <- {
  dd <- d[d$group == "All poor LS", ]
  ggplot(dd, aes(toxicity, efficacy)) + zone +
    geom_line(data = frontier(dd), aes(toxicity, efficacy), colour = "#7f8c8d", linetype = "dotted", linewidth = 0.55) +
    geom_point(aes(fill = label), shape = 21, size = 3.6, stroke = 0.3, colour = "grey20") +
    ggrepel::geom_text_repel(aes(label = add, colour = label), size = 2.6, fontface = "bold",
                             max.overlaps = 20, box.padding = 0.45, min.segment.length = 0,
                             segment.size = 0.2, seed = 1) +
    scale_fill_manual(values = TARGET12_COLORS, guide = "none") +
    scale_colour_manual(values = TARGET12_COLORS, guide = "none") +
    coord_cartesian(xlim = c(-4, 88), ylim = c(0, 100)) +
    labs(x = "normal HSPC / Myeloid toxicity (% covered)", y = "% of state cells covered (efficacy)") +
    ggtitle("All poor-prognosis LS (pooled)") + theme_pub(9)
}

ls_lv <- names(PROG)
lsd <- d |> dplyr::filter(group %in% ls_lv) |> dplyr::mutate(group = factor(group, levels = ls_lv))
frc <- lsd |> dplyr::group_by(group) |> dplyr::group_modify(~ frontier(.x)) |> dplyr::ungroup()
pC <- ggplot(lsd, aes(toxicity, efficacy)) + zone +
  geom_line(data = frc, aes(toxicity, efficacy), colour = "#7f8c8d", linetype = "dotted", linewidth = 0.4) +
  geom_point(aes(fill = label), shape = 21, size = 2.5, stroke = 0.2, colour = "grey25") +
  facet_wrap(~ group, ncol = 5, labeller = as_labeller(PROG)) +
  scale_fill_manual(values = TARGET12_COLORS, name = "Target",
                    guide = guide_legend(nrow = 1, override.aes = list(size = 6))) +
  coord_cartesian(xlim = c(-4, 100), ylim = c(0, 100)) +
  labs(x = "normal HSPC / Myeloid toxicity (% covered)", y = "% of state cells covered") +
  ggtitle("Per-state therapeutic window (CD96 + partner; dotted = Pareto frontier)") +
  theme_pub(8) +
  theme(strip.text = element_text(face = "bold"), legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 11), legend.text = element_text(size = 10),
        legend.key.size = grid::unit(0.45, "cm"))

fig <- pA / pC + plot_layout(heights = c(1, 0.95)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 15))
savefig(fig, "supplementary/FigureS23_LS_window", 12, 9)
cat("Done -> figures/supplementary/FigureS23_LS_window.{pdf,png}\n")
