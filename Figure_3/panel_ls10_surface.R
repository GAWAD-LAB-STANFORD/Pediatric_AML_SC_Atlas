#!/usr/bin/env Rscript
# LS10 surface-target panel: candidate HSC-sparing surface markers (present on leukemic
# cells incl. LS10, absent on normal HSC) across cell groups. Dot size = % expressing,
# fill = mean log-norm expression. Data from LS10_surface_dotplot.csv. No fabricated values.
suppressMessages({library(ggplot2); library(scales)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
D2 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
d <- read.csv(file.path(D2, "LS10_surface_dotplot.csv"))
d$group <- factor(d$group, levels=rev(c("normal HSC","other leukemic","LS3 (ns)","LS39 (ns)","LS10 (poor)")))
d$gene  <- factor(d$gene, levels=c("CD96","CLEC12A","SUCNR1","CD68"))
g <- ggplot(d, aes(gene, group, size=pct, fill=mean_expr)) +
  geom_point(shape=21, colour="grey30", stroke=0.3) +
  scale_size(range=c(0.5,11), limits=c(0,60), name="% expressing") +
  scale_fill_seq(name="mean\nlog-norm") +
  labs(x=NULL, y=NULL) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(face="bold", angle=0),
        axis.text.y=element_text(size=10), panel.grid.major=element_line(colour="grey92"),
        plot.caption=element_text(size=8, hjust=0))
ggsave(file.path(OUT, "Figure_2_LS10_surface.pdf"), g, width=5.2, height=3.2)
ggsave(file.path(OUT, "Figure_2_LS10_surface.png"), g, width=5.2, height=3.2, dpi=200)
cat("wrote Figure_2_LS10_surface\n")
