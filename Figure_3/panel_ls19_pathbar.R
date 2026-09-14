#!/usr/bin/env Rscript
# Figure 3 (LS19) pathway panel: horizontal bar chart of ALL significantly enriched GO:BP pathways
# (GSEA, LS19 vs all other leukemic cells; 0 significant DOWN, so every bar is an LS19-enriched term).
# Bar length = -log10 FDR (significance), fill = NES (enrichment) on the shared SEQ gradient (dot-plot
# scheme). Replaces the earlier word cloud (which tokenised term names and read as redundant).
# Data = LS19_gsea_sig.csv (from panel_ls19_pathways.R). No fabricated values.
suppressMessages({library(ggplot2); library(dplyr)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_3_panels"
d <- read.csv(file.path(D3, "LS19_gsea_sig.csv"))
d$nlp <- -log10(d$padj)
d <- d %>% arrange(nlp)                      # ascending -> most significant ends up on top
d$term <- factor(d$term, levels = d$term)
g <- ggplot(d, aes(x = nlp, y = term, fill = NES)) +
  geom_col(width = 0.72, colour = "grey30", linewidth = 0.2) +
  scale_fill_seq(name = "NES") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(-log[10]~FDR), y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8.5),
        legend.key.width = unit(0.35, "cm"))
ggsave(file.path(OUT, "Figure_3_LS19_pathbar.pdf"), g, width = 6.6, height = 5.0)
ggsave(file.path(OUT, "Figure_3_LS19_pathbar.png"), g, width = 6.6, height = 5.0, dpi = 200)
cat("wrote Figure_3_LS19_pathbar:", nrow(d), "significant UP pathways; FDR range",
    sprintf("%.0e to %.0e", min(d$padj), max(d$padj)), "\n")
