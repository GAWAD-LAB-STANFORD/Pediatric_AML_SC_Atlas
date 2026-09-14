#!/usr/bin/env Rscript
# Figure 3E (LS19) data-driven surface-target panel: the surfaceome genes expressed in the highest
# fraction of LS19 cells while nearly ABSENT on normal HSPC (therapeutic window that spares normal
# stem/progenitors), across normal HSPC / normal pDC / other leukemic / LS19. Dot size = % expressing,
# fill = mean log-norm (shared SEQ / dot-plot scheme). Genes ordered by LS19 % expressing; CD96 (the
# lead target and top hit) is highlighted. Data = compute_ls19_surface_screen.py. No fabricated values.
suppressMessages({library(ggplot2); library(scales)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_3_panels"
d <- read.csv(file.path(D3, "LS19_surface_screen.csv"), check.names = FALSE)
ls19 <- d[d$group == "LS19 (poor)", c("gene", "pct")]
gene_ord <- ls19$gene[order(-ls19$pct)]                       # CD96 first (highest LS19 %)
d$gene  <- factor(d$gene, levels = gene_ord)
d$group <- factor(d$group, levels = c("normal HSPC", "normal pDC", "other leukemic", "LS19 (poor)"))
cd96_x <- which(gene_ord == "CD96")
g <- ggplot(d, aes(gene, group, size = pct, fill = mean)) +
  annotate("rect", xmin = cd96_x - 0.5, xmax = cd96_x + 0.5, ymin = 0.4, ymax = 4.6,
           fill = "#9B2226", alpha = 0.08) +
  geom_point(shape = 21, colour = "grey30", stroke = 0.3) +
  scale_size(range = c(0.3, 9), name = "% expressing", labels = percent, limits = c(0, NA)) +
  scale_fill_seq(name = "mean\nlog-norm") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9,
          face = ifelse(gene_ord == "CD96", "bold", "plain")),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey92"),
        legend.key.size = unit(0.4, "cm"))
ggsave(file.path(OUT, "Figure_3_LS19_surface.pdf"), g, width = 7.2, height = 3.1)
ggsave(file.path(OUT, "Figure_3_LS19_surface.png"), g, width = 7.2, height = 3.1, dpi = 200)
cat("wrote Figure_3_LS19_surface (data-driven); genes:", paste(gene_ord, collapse = ", "), "\n")
