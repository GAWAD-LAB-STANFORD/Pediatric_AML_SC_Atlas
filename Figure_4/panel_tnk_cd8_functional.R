#!/usr/bin/env Rscript
# Figure 4 (T/NK): functional state across the CD8 differentiation axis (naive -> memory -> GZMK -> GZMB).
# The adverse, clonally-expanded GZMB-CD8 is maximally cytotoxic and terminally differentiated (high
# CX3CR1/FGFBP2/GZMH/KLRG1) but has lost stem/memory reserve (TCF7/IL7R/CD28 low) and is non-proliferative
# (MKI67 ~0), and is NOT classically exhausted (PD1/TIM3 low) - a spent terminal-effector, not exhaustion.
# Data = tnk_cd8_functional.csv. No fabricated values.
suppressMessages({library(ggplot2); library(dplyr); library(scales)})
source("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3/manuscript_palette.R")
D  <-"/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT<-"/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
d <- read.csv(file.path(D,"tnk_cd8_functional.csv"), check.names=FALSE)
prog_ord <- c("cytotoxic","terminal/senescence","memory/stem reserve","exhaustion","prolif")
d$program <- factor(d$program, levels=prog_ord)
d$gene <- factor(d$gene, levels=unique(d$gene[order(d$program)]))
d$subset <- factor(d$subset, levels=rev(c("Naïve CD8 T","Memory CD8 T","GZMK CD8 T","GZMB CD8 T")))
g <- ggplot(d, aes(gene, subset, size=pct, fill=mean)) +
  geom_point(shape=21, colour="grey30", stroke=0.3) +
  facet_grid(~program, scales="free_x", space="free_x", switch="x") +
  scale_size(range=c(0.3,7), name="% expressing", labels=percent) +
  scale_fill_seq(name="mean\nlog-norm") +
  labs(x=NULL, y=NULL) +
  theme_minimal(base_size=10) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8, face="italic"),
        strip.text=element_text(face="bold", size=8), strip.placement="outside",
        panel.grid.major=element_line(colour="grey93"), panel.spacing=unit(0.15,"lines"),
        legend.key.size=unit(0.35,"cm"))
ggsave(file.path(OUT,"Figure_4_cd8_functional.png"), g, width=9.0, height=2.8, dpi=200)
cat("wrote Figure_4_cd8_functional\n")
