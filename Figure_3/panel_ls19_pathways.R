#!/usr/bin/env Rscript
# Figure 3 (LS19) pathways panel: GSEA (fgsea, GO:BP) of LS19 vs other leukemic states (leukemic
# cells), genes ranked by Wilcoxon score. Only significant terms (BH p.adjust<0.05) shown. LS19 is
# non-cycling, so unlike LS10 the enriched pathways are a real program (interferon / antigen
# presentation / progenitor), not cell cycle. Also writes the significant-term table for the word
# cloud. Data = LS19_vs_other_DE.csv. No fabricated values.
suppressPackageStartupMessages({library(fgsea); library(msigdbr); library(dplyr); library(ggplot2); library(stringr)})
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_3_panels"
de <- read.csv(file.path(D3,"LS19_vs_other_DE.csv"))
r <- de$scores; names(r) <- de$names; r <- sort(r, decreasing=TRUE)
BP <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP") %>% dplyr::select(gs_name, gene_symbol)
paths <- split(BP$gene_symbol, BP$gs_name)
set.seed(42); fg <- fgseaMultilevel(paths, r, minSize=10, maxSize=500, eps=0)
sig <- fg[!is.na(fg$padj) & fg$padj<0.05 & fg$NES>0,]
sig$term <- tolower(gsub("_"," ", sub("^GOBP_","", sig$pathway)))
sig <- sig[order(sig$padj),]
write.csv(sig[,c("pathway","term","NES","padj","size")], file.path(D3,"LS19_gsea_sig.csv"), row.names=FALSE)
cat(sprintf("LS19 GSEA: %d significant up terms (padj<0.05)\n", nrow(sig)))
d <- head(sig, 12)
d$lab <- factor(str_wrap(d$term, 34), levels=rev(str_wrap(d$term,34)))
p <- ggplot(d, aes(-log10(padj), lab)) + geom_col(width=0.72, fill="#5aae61") +
  labs(x=expression(-log[10]~adjusted~p), y=NULL) +
  theme_classic(base_size=10) + theme(axis.text.y=element_text(size=8.5))
ggsave(file.path(OUT,"Figure_3_LS19_pathways.pdf"), p, width=5.8, height=4.6)
ggsave(file.path(OUT,"Figure_3_LS19_pathways.png"), p, width=5.8, height=4.6, dpi=200)
cat("top terms:", paste(head(d$term,8), collapse=" | "), "\n")
