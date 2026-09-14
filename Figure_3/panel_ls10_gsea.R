#!/usr/bin/env Rscript
# Figure 3C (significant-pathways version): GSEA (rank-based) of LS10 vs the other proliferation-
# high states (LS3, LS20, LS39), AML-leukemic cells only, genes ranked by Wilcoxon score. Only
# significant gene sets (BH p.adjust<0.05) are shown. NOTE: all significant pathways are enriched
# UP in LS10 (S-phase / replication / chromosome segregation); NO pathway is significantly enriched
# in the comparator proliferation states, so the panel is one-directional -- the differentiation
# difference is real at the gene level (LYZ/S100A4/6/LGALS1/TYMP down) but not a significant pathway.
# Data = LS10_vs_prolif_DE.csv. No fabricated values.
suppressPackageStartupMessages({library(clusterProfiler); library(msigdbr); library(dplyr); library(ggplot2); library(stringr)})
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
de <- read.csv(file.path(D3,"LS10_vs_prolif_DE.csv"))
r <- de$scores; names(r) <- de$names; r <- sort(r, decreasing=TRUE)
BP <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP") %>% dplyr::select(gs_name, gene_symbol)
set.seed(1)
g <- GSEA(r, TERM2GENE=BP, pvalueCutoff=0.05, minGSSize=10, maxGSSize=500, eps=0, verbose=FALSE)
res <- as.data.frame(g) %>% filter(p.adjust<0.05)
d <- res %>% filter(NES>0) %>% arrange(p.adjust) %>% head(10) %>%
  mutate(term=str_wrap(tolower(gsub("_"," ", sub("^GOBP_","", Description))), 34), nl=-log10(p.adjust))
d$term <- factor(d$term, levels=rev(d$term))
p <- ggplot(d, aes(nl, term)) +
  geom_col(width=0.72, fill="#d1495b") +
  labs(x=expression(-log[10]~adjusted~p), y=NULL) +
  theme_classic(base_size=10) +
  theme(axis.text.y=element_text(size=8.5))
ggsave(file.path(OUT,"Figure_3_LS10_gsea.pdf"), p, width=5.8, height=4.6)
ggsave(file.path(OUT,"Figure_3_LS10_gsea.png"), p, width=5.8, height=4.6, dpi=200)
cat(sprintf("wrote Figure_3_LS10_gsea: %d sig terms (%d up, %d down)\n", nrow(res), sum(res$NES>0), sum(res$NES<0)))
