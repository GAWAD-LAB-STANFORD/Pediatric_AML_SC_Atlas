#!/usr/bin/env Rscript
# Figure 3 (LS19) defining panel: marker programs across normal HSPC, LS19, normal pDC, and other
# leukemic cells (dot size = % expressing, colour = mean log-norm expression). Shows LS19 is an
# immature progenitor (SPINK2/PRSS57/CD34 high, like HSPC) with early pDC/interferon priming
# (IRF7/IFI6/partial IL3RA) but NOT mature pDC (LILRA4/CLEC4C absent) and low myeloid differentiation.
# Data = compute_ls19_markers.py. No fabricated values.
suppressMessages({library(ggplot2); library(dplyr); library(scales)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_3_panels"
d <- read.csv(file.path(D3,"LS19_marker_dotplot.csv"))
prog_ord <- c("progenitor","interferon","pDC","differentiation")   # surface -> its own panel (matches LS10)
d <- d[d$program %in% prog_ord,]
d$program <- factor(d$program, levels=prog_ord)
d <- d %>% arrange(program) %>% distinct(gene, group, .keep_all=TRUE)   # each gene -> its first program
d$group <- factor(d$group, levels=c("normal HSPC","LS19 (poor)","normal pDC","other leukemic"))
gene_ord <- unique(d$gene[order(d$program)])
d$gene <- factor(d$gene, levels=gene_ord)
p <- ggplot(d, aes(gene, group)) +
  geom_point(aes(size=pct, fill=mean), shape=21, colour="grey30", stroke=0.3) +
  facet_grid(~program, scales="free_x", space="free_x", switch="x") +
  scale_size(range=c(0.5,7), name="% expressing", labels=scales::percent) +
  scale_fill_seq(name="mean\nlog-norm") +
  labs(x=NULL, y=NULL) +
  theme_minimal(base_size=10) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8),
        strip.text=element_text(face="bold", size=8.5), strip.placement="outside",
        panel.grid.major=element_line(colour="grey93"), panel.spacing=unit(0.15,"lines"),
        legend.key.size=unit(0.35,"cm"))
ggsave(file.path(OUT,"Figure_3_LS19_program.pdf"), p, width=9, height=3.2)
ggsave(file.path(OUT,"Figure_3_LS19_program.png"), p, width=9, height=3.2, dpi=200)
cat("wrote Figure_3_LS19_program\n")
