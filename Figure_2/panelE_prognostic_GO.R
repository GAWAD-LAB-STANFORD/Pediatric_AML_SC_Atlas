#!/usr/bin/env Rscript
# Figure 2E (v3, corrected 49 states): GO:BP over-representation of the up-regulated
# genes per whole-cohort EFS-prognostic state (clusterProfiler::enricher vs MSigDB
# C5 GO:BP), rendered as a term x state heatmap (favorable | poor). Top-3 terms/state,
# unioned, capped at 28. Single-patient states marked with *. No fabricated values.
suppressPackageStartupMessages({library(clusterProfiler); library(msigdbr); library(dplyr)
  library(tidyr); library(ggplot2); library(scales); library(stringr); library(patchwork)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_2"
source(file.path(CODE, "manuscript_palette.R"))
D2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
UP  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3/LS_upgenes_wilcoxon.csv"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"

get_bp <- function() for (call in list(
    function() msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP"),
    function() msigdbr(species="Homo sapiens", category="C5", subcategory="GO:BP"),
    function() msigdbr(species="Homo sapiens", category="C5", subcategory="BP")))
  { r <- tryCatch(call(), error=function(e) NULL); if (!is.null(r) && nrow(r)) return(r) }
BP <- get_bp() %>% dplyr::select(gs_name, gene_symbol)

# order: favorable first, then poor, each by p; labels with program + single-patient *
pr <- read.csv(file.path(D2, "LS_wholecohort_prognosis.csv")); pr <- pr[pr$prognosis %in% c("favorable","poor"), ]
pr$dir <- pr$prognosis; pr$star <- ifelse(pr$n_patients <= 1, "*", "")
pr$progn <- ifelse(pr$lineage_program=="unresolved", "", pr$lineage_program)
pr$lab <- ifelse(pr$progn=="", paste0(gsub("_","",pr$LS),pr$star), paste0(gsub("_","",pr$LS),pr$star,"  ",pr$progn))
ord <- pr[order(pr$dir=="poor", pr$p), ]; ord$lab <- make.unique(ord$lab)

up <- read.csv(UP)
# NORMALIZE FOR GENE NUMBER: give every state an equal-size input gene set = its top-N
# up-genes by log2FC, N = the smallest prognostic state's up-gene count. Removes the
# over-representation gene-set-size confound (bigger up-lists hit more/stronger terms).
TOPN <- min(table(factor(up$LS[up$LS %in% ord$LS], levels=ord$LS)))
cat(sprintf("equal-N input: top-%d up-genes (by log2FC) per state\n", TOPN))
res <- list()
for (i in seq_len(nrow(ord))) {
  x <- ord$LS[i]; ug <- up[up$LS == x, ]; ug <- ug[order(-ug$log2FC), ]
  g <- head(ug$gene, TOPN); if (length(g) < 5) next
  em <- tryCatch(enricher(g, TERM2GENE=BP, pvalueCutoff=1, qvalueCutoff=1, minGSSize=5, maxGSSize=500), error=function(e) NULL)
  if (is.null(em) || !nrow(as.data.frame(em))) next
  d <- as.data.frame(em); d$LS <- x; res[[x]] <- d[, c("LS","Description","p.adjust")]
}
R <- bind_rows(res)
R$term <- R$Description |> sub("^GOBP_","",x=_) |> gsub("_"," ",x=_) |> tolower()
top <- R %>% group_by(LS) %>% slice_min(p.adjust, n=3, with_ties=FALSE) %>% ungroup()
terms <- unique(top$term)
if (length(terms) > 28) terms <- R %>% filter(term %in% terms) %>% group_by(term) %>%
  summarise(m=min(p.adjust)) %>% slice_min(m, n=28) %>% pull(term)
M <- R %>% filter(term %in% terms) %>% mutate(nl=-log10(p.adjust))
full <- expand.grid(term=terms, LS=ord$LS, stringsAsFactors=FALSE) %>%
  left_join(M[,c("LS","term","nl")], by=c("LS","term")) %>% mutate(nl=ifelse(is.na(nl),0,nl)) %>%
  left_join(ord[,c("LS","lab","dir")], by="LS")
peak <- full %>% group_by(term) %>% slice_max(nl, n=1, with_ties=FALSE) %>%
  mutate(ci=match(LS, ord$LS)) %>% select(term, ci) %>% arrange(ci)
full$term <- factor(full$term, levels=rev(peak$term))
full$lab  <- factor(full$lab,  levels=ord$lab)
full$term <- factor(full$term, levels=levels(full$term), labels=str_wrap(levels(full$term), 34))
nfav <- sum(ord$dir=="favorable"); xdiv <- nfav + 0.5
labcol <- ifelse(ord$dir=="favorable", unname(PROG["favorable"]), unname(PROG["poor"]))

hp <- ggplot(full, aes(lab, term, fill=nl)) +
  geom_tile(colour="white", linewidth=0.4) +
  geom_vline(xintercept=xdiv, colour="grey30", linewidth=0.6) +
  scale_fill_seq(limits=c(0,6), oob=squish, name=expression(-log[10]~q)) +
  labs(x=NULL, y=NULL) + theme_minimal(base_size=9) +
  theme(axis.text.x=element_text(angle=45, hjust=1, face="bold", size=8.5, colour=labcol),
        axis.text.y=element_text(size=7.5), panel.grid=element_blank(),
        plot.caption=element_text(size=7.5, hjust=0), plot.margin=margin(0,10,6,6))
# stemness strip (HSC/MPP program z, same metric as the 2C strip), aligned above the heatmap
strip <- ord; strip$lab <- factor(strip$lab, levels=ord$lab)
sp <- ggplot(strip, aes(lab, "stemness", fill=stem_z)) +
  geom_tile(colour="white", linewidth=0.4) +
  geom_vline(xintercept=xdiv, colour="grey30", linewidth=0.6) +
  scale_fill_div(midpoint=0, limits=c(-1,1), oob=squish, name="stemness (z)") +
  labs(x=NULL, y=NULL) + theme_minimal(base_size=9) +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=8),
        panel.grid=element_blank(), plot.margin=margin(6,10,0,6))
fig <- sp / hp + plot_layout(heights=c(0.045,1), guides="collect")
wd <- max(6.5, 3.4 + 0.42*nrow(ord))
ggsave(file.path(OUT, "Figure_2E_prognostic_GOBP.pdf"), fig, width=wd, height=7.6)
ggsave(file.path(OUT, "Figure_2E_prognostic_GOBP.png"), fig, width=wd, height=7.6, dpi=200)
write.csv(R, file.path(D2, "LS_prognostic_GOBP_all.csv"), row.names=FALSE)
cat(sprintf("wrote Figure_2E_prognostic_GOBP: %d terms x %d states\n", length(terms), nrow(ord)))
