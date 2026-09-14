#!/usr/bin/env Rscript
# Figure 3C: the 162 GO:BP pathways significantly enriched in LS10 (GSEA vs the other proliferation-
# high states LS3/LS20/LS39, leukemic cells only; BH p<0.05, all up in LS10) collapsed into logical
# themes. Each theme is assigned by priority keyword; the bar is the strongest term in the theme
# (-log10 adjusted p), and the number of significant terms in the theme is annotated. All themes are
# replication/mitosis/repair (the significant signal is the replication program); the differentiation
# difference is gene-level only. Data = LS10_vs_prolif_DE.csv. No fabricated values.
suppressPackageStartupMessages({library(fgsea); library(msigdbr); library(dplyr); library(ggplot2)})
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
de <- read.csv(file.path(D3,"LS10_vs_prolif_DE.csv"))
BP <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP"); paths <- split(BP$gene_symbol, BP$gs_name)
s <- de$scores; names(s) <- de$names; s <- sort(s, decreasing=TRUE)
set.seed(42); fg <- fgseaMultilevel(paths, s, minSize=10, maxSize=500, eps=0)
sig <- fg[fg$padj<0.05 & fg$NES>0,]
sig$term <- tolower(gsub("_"," ", sub("^GOBP_","", sig$pathway)))
THEMES <- list(
  "DNA replication"            = "replicat",
  "Chromosome segregation"     = "segregation|chromatid|kinetochore|centromere|chromosome separation|condensation",
  "Mitosis & nuclear division" = "mitotic|mitosis|nuclear division|spindle|cytokinesis|metaphase|microtubule",
  "Cell-cycle regulation"      = "cell cycle|phase transition|checkpoint|g1|g2|s phase|g0",
  "DNA repair & recombination" = "repair|recombinat|damage",
  "Chromatin & telomere"       = "chromatin|histone|nucleosome|telomere")
assign_theme <- function(x){ for(nm in names(THEMES)) if(grepl(THEMES[[nm]], x)) return(nm); "Other" }
sig$theme <- factor(sapply(sig$term, assign_theme), levels=c(names(THEMES),"Other"))
grp <- sig %>% group_by(theme) %>% summarise(n=n(), nl=max(-log10(padj)), .groups="drop") %>%
  filter(!is.na(theme)) %>% arrange(nl)
grp$theme <- factor(grp$theme, levels=grp$theme)
write.csv(grp, file.path(D3,"LS10_gsea_themes.csv"), row.names=FALSE)
p <- ggplot(grp, aes(nl, theme)) +
  geom_col(width=0.7, fill="#d1495b") +
  geom_text(aes(label=paste0(n," terms")), hjust=-0.15, size=3, colour="grey30") +
  scale_x_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(x=expression(-log[10]~adjusted~p~(strongest~term)), y=NULL) +
  theme_classic(base_size=11) + theme(axis.text.y=element_text(size=9.5))
ggsave(file.path(OUT,"Figure_3_LS10_gsea.pdf"), p, width=5.8, height=3.6)
ggsave(file.path(OUT,"Figure_3_LS10_gsea.png"), p, width=5.8, height=3.6, dpi=200)
cat("themes:\n"); print(as.data.frame(grp))
PY
