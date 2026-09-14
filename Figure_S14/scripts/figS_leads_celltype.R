# ====================================================================
# figS_leads_celltype.R  |  Figure S16 (lead-target cell-type specificity)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : figS16_leads_celltype.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript figS_leads_celltype.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# ALTERNATE — cell-type specificity / safety of the 4 LEAD targets (CD96, CD9, SUCNR1, IL1RAP).
# % of cells positive across normal compartments vs core leukemic. The funnel's toxicity gate only
# inspects 0_HSPC + Myeloid_Pro, so this is the fuller safety view. Source: figS16_leads_celltype.csv.
source("_style.R")

d <- sd("figS16_leads_celltype.csv")
comp_lv <- c("leukemic", "HSPC / Myeloid (gate)", "monocyte / DC", "T / NK", "B / plasma", "erythroid")
d$compartment <- factor(d$compartment, levels = comp_lv)
d$gene <- factor(d$gene, levels = LEAD_TARGETS)
# cell-type order within each compartment (by total expression, leukemic kept on top)
ord <- d |> dplyr::group_by(cell_type) |> dplyr::summarise(s = sum(pct)) |> dplyr::arrange(s)
d$cell_type <- factor(d$cell_type, levels = ord$cell_type)

gx <- function(v) sprintf("<b style='color:%s'>%s</b>", LEAD_COLORS[as.character(v)], v)

p <- ggplot(d, aes(gene, cell_type)) +
  geom_point(aes(size = pct, fill = pct), shape = 21, colour = "#3b5870", stroke = 0.25) +
  scale_fill_gradientn(colours = SEQ, limits = c(0, 65), oob = scales::squish, name = "% positive") +
  scale_size_area(max_size = 8, limits = c(0, 65), name = "% positive", breaks = c(5, 20, 40, 60)) +
  scale_x_discrete(labels = function(v) gx(v)) +
  facet_grid(compartment ~ ., scales = "free_y", space = "free", switch = "y") +
  labs(x = NULL, y = NULL,
       title = "Cell-type specificity of the four lead targets",
       subtitle = NULL) +
  theme_pub(9) +
  theme(axis.text.x = ggtext::element_markdown(size = 11),
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7.5, face = "bold"),
        strip.placement = "outside", panel.spacing = unit(2, "pt"),
        plot.subtitle = ggtext::element_textbox_simple(size = 8, colour = "#34495e", margin = margin(b = 6)))

savefig(p, "supplementary/FigureS16_leads_celltype_specificity", 6.5, 8)
cat("Done -> figures/supplementary/FigureS16_leads_celltype_specificity.{pdf,png}\n")
