#!/usr/bin/env Rscript
# Figure 4 pathway heatmap (per-subset GSEA, then combine): union of the 5 most significant GO:BP pathways
# per T/NK subset, coloured by signed -log10 FDR (sign of NES; red = up in AML, blue = up in normal) on the
# manuscript diverging palette. Cells where fgsea could not estimate a reliable FDR for a subset are filled
# as neutral 0 (not significant) rather than left blank. Data = tnk_gsea_persubset.csv. No fabrication.
suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})
source("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3/manuscript_palette.R")
V2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
ORD <- c("Naïve CD4 T","Memory CD4 T","Treg","Naïve CD8 T","Memory CD8 T","GZMK CD8 T",
         "GZMB CD8 T","GZMB DNT","MAIT","GZMK NK","GZMB NK","Proliferating T")
all <- read.csv(file.path(V2, "tnk_gsea_persubset.csv"), check.names = FALSE)
top5 <- all %>% filter(padj < 0.05) %>% group_by(subset) %>% slice_min(padj, n = 5, with_ties = FALSE) %>% ungroup()
keep <- unique(top5$term)
mat <- all %>% filter(term %in% keep) %>% mutate(signif = sign(NES) * pmin(-log10(padj), 10)) %>%
  select(subset, term, signif) %>%
  complete(subset = ORD, term = keep, fill = list(signif = 0))      # fill un-estimated cells as neutral 0
mat$subset <- factor(mat$subset, levels = ORD)
ord_p <- top5 %>% count(term) %>% arrange(desc(n)) %>% pull(term)
ord_p <- unique(c(ord_p, setdiff(keep, ord_p)))
mat$term <- factor(mat$term, levels = rev(ord_p))
g <- ggplot(mat, aes(subset, term, fill = signif)) +
  geom_tile(colour = "grey88") +
  scale_fill_div(midpoint = 0, name = "signed\n-log10 FDR", limits = c(-10, 10)) +
  labs(x = NULL, y = NULL, caption = "red = enriched in AML, blue = enriched in normal; union of top-5 pathways per subset") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 7), panel.grid = element_blank(),
        plot.caption = element_text(size = 8, hjust = 0))
ggsave(file.path(OUT, "Figure_4H_tnk_gsea_heatmap.pdf"), g, width = 8.0, height = 9.0)
ggsave(file.path(OUT, "Figure_4H_tnk_gsea_heatmap.png"), g, width = 8.0, height = 9.0, dpi = 200)
cat("wrote Figure_4H_tnk_gsea_heatmap (", length(keep), "pathways, grid completed )\n")
