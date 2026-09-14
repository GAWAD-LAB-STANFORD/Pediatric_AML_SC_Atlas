# ====================================================================
# figS_cyto_window.R  |  Figure S24 (additive window by cytogenetic subtype)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : figS_cyto_window.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript figS_cyto_window.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# Supplementary figure — therapeutic window by CYTOGENETIC SUBTYPE (analogue of S15, grouped
# by cytogenetics instead of PAC). Each point = CD96 alone or CD96 OR one partner; x = normal
# HSPC/Myeloid toxicity, y = % of that subtype's leukemic cells covered; points sorted by
# toxicity, Pareto frontier (dotted). Source: figS_cyto_window.csv (build_cyto_window.py).
source("_style.R")

d <- sd("figS_cyto_window.csv")
d$label <- factor(d$label, levels = names(TARGET12_COLORS))
GREEN <- "#eafaf1"; RED <- "#fdecea"
zone <- list(annotate("rect", xmin = -4, xmax = 10, ymin = -Inf, ymax = Inf, fill = GREEN),
             annotate("rect", xmin = 10, xmax = 115, ymin = -Inf, ymax = Inf, fill = RED))
frontier <- function(dd) { dd <- dd[order(dd$toxicity), ]; dd[dd$efficacy >= cummax(dd$efficacy) - 1e-9, ] }

## A — All AML (pooled leukemic), labeled
dd <- d |> dplyr::filter(group == "All AML")
pA <- ggplot(dd, aes(toxicity, efficacy)) + zone +
  geom_line(data = frontier(dd), aes(toxicity, efficacy), colour = "#7f8c8d", linetype = "dotted", linewidth = 0.55) +
  geom_point(aes(fill = label), shape = 21, size = 3.6, stroke = 0.3, colour = "grey20") +
  ggrepel::geom_text_repel(aes(label = add, colour = label), size = 2.5, fontface = "bold",
                           max.overlaps = 20, box.padding = 0.45, min.segment.length = 0, seed = 1) +
  scale_fill_manual(values = TARGET12_COLORS, guide = "none") +
  scale_colour_manual(values = TARGET12_COLORS, guide = "none") +
  coord_cartesian(xlim = c(-4, 88), ylim = c(0, 100)) +
  labs(x = "normal HSPC / Myeloid toxicity (% covered)", y = "% of leukemic cells covered (efficacy)",
       title = "All AML (pooled leukemic cells)") +
  theme_pub(8.5)

## B — per-cytogenetic-subtype facets (ordered by CD96-alone efficacy) + big key
cyto_ord <- d |> dplyr::filter(group_type == "cyto", add == "CD96 alone") |>
  dplyr::arrange(dplyr::desc(efficacy)) |> dplyr::pull(group)
cytod <- d |> dplyr::filter(group_type == "cyto") |> dplyr::mutate(group = factor(group, levels = cyto_ord))
frc <- cytod |> dplyr::group_by(group) |> dplyr::group_modify(~ frontier(.x)) |> dplyr::ungroup()
pB <- ggplot(cytod, aes(toxicity, efficacy)) + zone +
  geom_line(data = frc, aes(toxicity, efficacy), colour = "#7f8c8d", linetype = "dotted", linewidth = 0.4) +
  geom_point(aes(fill = label), shape = 21, size = 2.3, stroke = 0.2, colour = "grey25") +
  facet_wrap(~ group, ncol = 5) +
  scale_fill_manual(values = TARGET12_COLORS, name = "Target",
                    guide = guide_legend(nrow = 1, override.aes = list(size = 6))) +
  coord_cartesian(xlim = c(-4, 100), ylim = c(0, 100)) +
  labs(x = "normal HSPC / Myeloid toxicity (% covered)", y = "% of leukemic cells covered",
       title = "Per-cytogenetic-subtype therapeutic window (CD96 + partner; dotted = frontier)") +
  theme_pub(8) +
  theme(strip.text = element_text(face = "bold", size = 6.8), legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11), legend.key.size = grid::unit(0.5, "cm"))

fig <- pA / pB + plot_layout(heights = c(1, 2.4)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))
savefig(fig, "supplementary/FigureS24_cyto_window", 12, 15)
