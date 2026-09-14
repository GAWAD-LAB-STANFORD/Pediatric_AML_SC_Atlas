#!/usr/bin/env Rscript
# LS10 pathway panel: GO:BP over-representation of genes UP in LS10 vs UP in the NS-
# proliferation comparators (LS39, LS3). Two facets, -log10(padj). Shows LS10 = intense
# undifferentiated proliferation (DNA replication / chromosome segregation), comparators =
# myeloid/immune differentiation. Gene lists from the per-cell DE. No fabricated values.
suppressPackageStartupMessages({library(clusterProfiler); library(msigdbr); library(dplyr); library(ggplot2); library(stringr)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
D2 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
BP <- msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP") %>% dplyr::select(gs_name, gene_symbol)
goterms <- function(csv, n=7){
  g <- read.csv(file.path(D2, csv))$names
  em <- as.data.frame(enricher(g, TERM2GENE=BP, pvalueCutoff=1, qvalueCutoff=1, minGSSize=5, maxGSSize=500))
  em$term <- tolower(gsub("_"," ", sub("^GOBP_","", em$Description)))
  em %>% arrange(p.adjust) %>% head(n) %>% transmute(term, nl=-log10(p.adjust))
}
up  <- goterms("LS10_up_nonhistone.csv"); up$dir  <- "up in LS10 (proliferation)"
dn  <- goterms("LS10_dn.csv");            dn$dir  <- "up in comparator (differentiation)"
d <- bind_rows(up, dn)
d$dir <- factor(d$dir, levels=c("up in LS10 (proliferation)","up in comparator (differentiation)"))
d <- d %>% group_by(dir) %>% mutate(term=str_wrap(term, 28)) %>% ungroup()
d$term <- factor(d$term, levels=rev(d$term))
g <- ggplot(d, aes(nl, term, fill=dir)) +
  geom_col(width=0.72) +
  facet_wrap(~dir, scales="free", ncol=1) +
  scale_fill_manual(values=c("up in LS10 (proliferation)"=unname(PROG["poor"]),
                             "up in comparator (differentiation)"=unname(PROG["favorable"])), guide="none") +
  labs(x=expression(-log[10]~adjusted~p), y=NULL) +
  theme_classic(base_size=10) +
  theme(strip.text=element_text(face="bold", size=9.5), strip.background=element_blank(),
        axis.text.y=element_text(size=8.5))
ggsave(file.path(OUT, "Figure_2_LS10_pathways.pdf"), g, width=5.4, height=5.0)
ggsave(file.path(OUT, "Figure_2_LS10_pathways.png"), g, width=5.4, height=5.0, dpi=200)
cat("wrote Figure_2_LS10_pathways\n")
