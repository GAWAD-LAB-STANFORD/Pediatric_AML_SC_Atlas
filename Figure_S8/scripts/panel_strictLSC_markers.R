#!/usr/bin/env Rscript
# Panel: the 26 strict-LSC candidates are NOT primitive quiescent LSC -- they co-express
# GMP/granule differentiation markers (AZU1/CTSG/MPO), placing them between primitive HSC
# and committed GMP. Compared to normal HSC (primitive) and normal GMP (committed) across
# four marker modules. Dot size = % expressing, fill = mean log-CPM. No fabricated values.
suppressMessages({library(ggplot2); library(scales)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
D3 <- "../source_data"
OUT <- "../figures/Figure_2_panels"
d <- read.csv(file.path(D3, "strictLSC_marker_profile.csv"))
d$group  <- factor(d$group, levels = rev(c("normal HSC (primitive)","strict-LSC candidates (n=26)","normal GMP (committed)")))
d$module <- factor(d$module, levels = c("stem/quiescence","LSC-aberrant","GMP/granule","monocyte"))
d$gene   <- factor(d$gene, levels = c("CD34","HLF","AVP","HOXA9","CD96","CLEC12A","HAVCR2","MPO","ELANE","AZU1","CTSG","LYZ","CD14"))
g <- ggplot(d, aes(gene, group, size = pct, fill = mean_expr)) +
  geom_point(shape = 21, colour = "grey30", stroke = 0.3) +
  facet_grid(~module, scales = "free_x", space = "free_x") +
  scale_size(range = c(0.3, 9), limits = c(0, 100), name = "% expressing") +
  scale_fill_seq(name = "mean\nlog-CPM") +
  labs(x = NULL, y = NULL,
       caption = "The strict-LSC candidates carry GMP/granule differentiation markers (AZU1/CTSG/MPO) -- aberrant GMP-like cells, not primitive quiescent LSC.") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(face = "bold", size = 8.5), axis.text.y = element_text(size = 9.5),
        strip.text = element_text(face = "bold", size = 9), panel.grid.major = element_line(colour = "grey92"),
        plot.caption = element_text(size = 7.4, hjust = 0))
ggsave(file.path(OUT, "Figure_3_strictLSC_markers.pdf"), g, width = 8.6, height = 3.2)
ggsave(file.path(OUT, "Figure_3_strictLSC_markers.png"), g, width = 8.6, height = 3.2, dpi = 200)
cat("wrote Figure_3_strictLSC_markers\n")
