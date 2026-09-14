#!/usr/bin/env Rscript
# Figure 4: GSEA of the AML-vs-healthy-BM difference in T/NK cells. Genes are ranked by the mean
# within-subset logFC across the 12 T/NK subsets (composition-controlled AML-vs-HBM axis). fgsea on GO:BP.
# HONEST NOTE: the ranking is dominated by an immediate-early/stress (dissociation) artifact, so we also
# re-run after removing the canonical artifact gene classes to see what real biology survives. No fabrication.
suppressMessages({library(fgsea); library(msigdbr); library(dplyr)})
V2 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
de <- read.csv(file.path(V2, "tnk_aml_vs_hbm_de_bySubset.csv"), check.names = FALSE)
rk <- de %>% group_by(gene) %>% summarise(lfc = mean(logFC), .groups = "drop")
BP <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
      dplyr::select(gs_name, gene_symbol)
paths <- split(BP$gene_symbol, BP$gs_name)
run_gsea <- function(rkdf, tag) {
  r <- rkdf$lfc; names(r) <- rkdf$gene; r <- sort(r, decreasing = TRUE)
  set.seed(42)
  fg <- as.data.frame(fgseaMultilevel(paths, r, minSize = 10, maxSize = 500, eps = 0))
  fg$term <- tolower(gsub("_", " ", sub("^GOBP_", "", fg$pathway)))
  fg <- fg[!is.na(fg$padj), ]
  cat(sprintf("\n=== %s: %d sig pathways (padj<0.05): UP(AML) %d / DOWN(normal) %d ===\n",
      tag, sum(fg$padj < 0.05), sum(fg$padj < 0.05 & fg$NES > 0), sum(fg$padj < 0.05 & fg$NES < 0)))
  up <- fg %>% filter(padj < 0.05, NES > 0) %>% arrange(padj) %>% head(15)
  dn <- fg %>% filter(padj < 0.05, NES < 0) %>% arrange(padj) %>% head(15)
  cat("-- UP in AML (top 15) --\n"); print(up[, c("term","NES","padj")], row.names = FALSE)
  cat("-- DOWN (up in normal, top 15) --\n"); print(dn[, c("term","NES","padj")], row.names = FALSE)
  fg
}
fg_all <- run_gsea(rk, "ALL genes")
write.csv(fg_all[order(fg_all$padj), c("term","NES","padj","size")], file.path(V2, "tnk_gsea_all.csv"), row.names = FALSE)
# de-artifacted: drop IEG/stress/heat-shock/ribosomal/Ig classes
IEG <- c("JUN","JUNB","JUND","FOS","FOSB","FOSL1","FOSL2","EGR1","EGR2","EGR3","IER2","IER3","IER5","ATF3",
 "NR4A1","NR4A2","NR4A3","DUSP1","DUSP2","ZFP36","ZFP36L1","ZFP36L2","BTG2","GADD45B","PPP1R15A","KLF2","KLF4","KLF6",
 "YPEL5","SOCS3","NFKBIA","NFKBIZ","CD69","MCL1","RGS1","CEBPD","CEBPB","IRF1","TNFAIP3","UBC")
art <- function(g) g %in% IEG | grepl("^(HSP|RPS|RPL|MT-|MTRNR|IG[KHL][VC]|DNAJ)", g)
fg_bio <- run_gsea(rk %>% filter(!art(gene)), "ARTIFACT-REMOVED")
write.csv(fg_bio[order(fg_bio$padj), c("term","NES","padj","size")], file.path(V2, "tnk_gsea_deartifacted.csv"), row.names = FALSE)
