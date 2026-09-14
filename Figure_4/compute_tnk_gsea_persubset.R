#!/usr/bin/env Rscript
# Figure 4: per-subset GSEA then combine (the rigorous version). fgsea GO:BP is run SEPARATELY within each
# T/NK subset on that subset's AML-vs-HBM logFC ranking, so we can see whether pathways are shared across
# subsets or driven by one. Heatmap = the union of the 5 most significant pathways per subset, coloured by
# signed significance (sign of NES x -log10 padj: red = up in AML, blue = up in normal). Data =
# tnk_aml_vs_hbm_de_bySubset.csv. No fabricated values.
suppressMessages({library(fgsea); library(msigdbr); library(dplyr); library(ggplot2); library(tidyr)})
source("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3/manuscript_palette.R")
V2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
de <- read.csv(file.path(V2, "tnk_aml_vs_hbm_de_bySubset.csv"), check.names = FALSE)
BP <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP") %>% dplyr::select(gs_name, gene_symbol)
paths <- split(BP$gene_symbol, BP$gs_name)
ORD <- c("Naïve CD4 T","Memory CD4 T","Treg","Naïve CD8 T","Memory CD8 T","GZMK CD8 T",
         "GZMB CD8 T","GZMB DNT","MAIT","GZMK NK","GZMB NK","Proliferating T")
subs <- intersect(ORD, unique(de$subset))
res <- list()
for (s in subs) {
  d <- de[de$subset == s, ]; r <- d$logFC; names(r) <- d$gene; r <- sort(r, decreasing=TRUE)
  set.seed(42)
  fg <- as.data.frame(fgseaMultilevel(paths, r, minSize=10, maxSize=500, eps=0))
  fg$subset <- s; fg$term <- tolower(gsub("_"," ", sub("^GOBP_","", fg$pathway)))
  res[[s]] <- fg[!is.na(fg$padj), ]
  cat(sprintf("%-16s sig pathways padj<0.05: %d (up %d / down %d)\n", s, sum(fg$padj<0.05,na.rm=TRUE),
      sum(fg$padj<0.05 & fg$NES>0,na.rm=TRUE), sum(fg$padj<0.05 & fg$NES<0,na.rm=TRUE)))
}
all <- bind_rows(res)
write.csv(all[,c("subset","term","NES","padj","size")], file.path(V2,"tnk_gsea_persubset.csv"), row.names=FALSE)
# top 5 significant pathways per subset -> union
top5 <- all %>% filter(padj<0.05) %>% group_by(subset) %>% slice_min(padj, n=5, with_ties=FALSE) %>% ungroup()
keep <- unique(top5$term)
mat <- all %>% filter(term %in% keep) %>%
  mutate(signif = sign(NES) * pmin(-log10(padj), 10)) %>%          # signed, capped at 10 for display
  select(subset, term, signif)
mat$subset <- factor(mat$subset, levels = ORD)
# order pathways by how many subsets they're significant in (then mean signif)
ord_p <- top5 %>% count(term) %>% arrange(desc(n)) %>% pull(term)
ord_p <- unique(c(ord_p, setdiff(keep, ord_p)))
mat$term <- factor(mat$term, levels = rev(ord_p))
g <- ggplot(mat, aes(subset, term, fill = signif)) +
  geom_tile(colour = "grey88") +
  scale_fill_div(midpoint = 0, name = "signed\n-log10 FDR", limits = c(-10,10)) +
  labs(x = NULL, y = NULL, caption = "red = enriched in AML, blue = enriched in normal; union of top-5 pathways per subset") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 7), panel.grid = element_blank(),
        plot.caption = element_text(size = 8, hjust = 0))
ggsave(file.path(OUT, "Figure_4H_tnk_gsea_heatmap.pdf"), g, width = 8.0, height = 9.0)
ggsave(file.path(OUT, "Figure_4H_tnk_gsea_heatmap.png"), g, width = 8.0, height = 9.0, dpi = 200)
cat("\nunion of top-5-per-subset pathways:", length(keep), "\nwrote Figure_4H_tnk_gsea_heatmap + tnk_gsea_persubset.csv\n")
