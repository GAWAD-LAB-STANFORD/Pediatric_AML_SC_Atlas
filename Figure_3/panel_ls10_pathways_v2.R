#!/usr/bin/env Rscript
# Figure 3C: GO:BP over-representation for genes distinguishing LS10 from the OTHER proliferation-
# high states (LS3, LS20, LS39) -- this cancels proliferation itself. Two facets: pathways up in
# LS10 (E2F/S-phase replication; histone genes excluded so terms reflect regulators not histones)
# and pathways up in the other proliferation states (myeloid differentiation / immune). Message:
# among proliferating blasts LS10 uniquely couples replication to a differentiation block.
# Gene lists from LS10_vs_prolif_DE.csv (Wilcoxon, BH-adjusted). No fabricated values.
suppressPackageStartupMessages({library(clusterProfiler); library(msigdbr); library(dplyr); library(ggplot2); library(stringr)})
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
de <- read.csv(file.path(D3,"LS10_vs_prolif_DE.csv"))
up <- de %>% filter(pvals_adj<0.05, logfoldchanges> 0.5, !grepl("^HIST", names)) %>% pull(names)
dn <- de %>% filter(pvals_adj<0.05, logfoldchanges< -0.5) %>% pull(names)
BP <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP") %>% dplyr::select(gs_name, gene_symbol)
goterms <- function(g, n=7){
  em <- as.data.frame(enricher(g, TERM2GENE=BP, pvalueCutoff=1, qvalueCutoff=1, minGSSize=5, maxGSSize=500))
  if(nrow(em)==0) return(NULL)
  em$term <- tolower(gsub("_"," ", sub("^GOBP_","", em$Description)))
  em %>% arrange(p.adjust) %>% head(n) %>% transmute(term, nl=-log10(p.adjust))
}
u <- goterms(up); u$dir <- "up in LS10 (replication)"
d <- goterms(dn); d$dir <- "up in other proliferation states (differentiation)"
dd <- bind_rows(u,d)
dd$dir <- factor(dd$dir, levels=c("up in LS10 (replication)","up in other proliferation states (differentiation)"))
dd <- dd %>% group_by(dir) %>% mutate(term=str_wrap(term,30)) %>% ungroup()
dd$term <- factor(dd$term, levels=rev(unique(dd$term)))
g <- ggplot(dd, aes(nl, term, fill=dir)) +
  geom_col(width=0.72) + facet_wrap(~dir, scales="free", ncol=1) +
  scale_fill_manual(values=c("up in LS10 (replication)"="#d1495b",
                             "up in other proliferation states (differentiation)"="#00798c"), guide="none") +
  labs(x=expression(-log[10]~adjusted~p), y=NULL,
       title="Pathways distinguishing LS10 among proliferating states") +
  theme_classic(base_size=10) +
  theme(strip.text=element_text(face="bold", size=9), strip.background=element_blank(),
        axis.text.y=element_text(size=8), plot.title=element_text(size=10.5,face="bold"))
ggsave(file.path(OUT,"Figure_3_LS10_pathways.pdf"), g, width=5.6, height=5.2)
ggsave(file.path(OUT,"Figure_3_LS10_pathways.png"), g, width=5.6, height=5.2, dpi=200)
cat("wrote Figure_3_LS10_pathways: up terms=",nrow(u)," down terms=",nrow(d),"\n")
