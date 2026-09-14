#!/usr/bin/env Rscript
# Panel G (GO:BP version): GO biological-process over-representation of the
# up-regulated genes per prognosis-associated state. Method mirrors the original
# Figure 2F (clusterProfiler::enricher vs MSigDB C5 GO:BP), rendered as a heatmap.
suppressPackageStartupMessages({library(clusterProfiler); library(msigdbr); library(dplyr)
  library(tidyr); library(ggplot2); library(patchwork); library(scales); library(stringr)})
SL <- "/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
source(file.path(SL, "manuscript_palette.R"))
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_new_panels"
args <- commandArgs(trailingOnly = TRUE)
ORDER <- if (length(args) >= 1) args[1] else file.path(SL, "ls48/LS_prognostic_order.csv")
OUTB  <- if (length(args) >= 2) args[2] else "Figure_2G_prognostic_GOBP"

get_bp <- function() for (call in list(
    function() msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP"),
    function() msigdbr(species="Homo sapiens", category="C5", subcategory="GO:BP"),
    function() msigdbr(species="Homo sapiens", category="C5", subcategory="BP")))
  { r <- tryCatch(call(), error=function(e) NULL); if (!is.null(r) && nrow(r)) return(r) }
BP <- get_bp() %>% dplyr::select(gs_name, gene_symbol)

up  <- read.csv(file.path(SL, "ls48/LS48_upgenes_wilcoxon.csv"))
ord <- read.csv(ORDER)                                               # LS, lab, group, dir, p (fav-first)
res <- list()
for (i in seq_len(nrow(ord))) {
  x <- ord$LS[i]; g <- up$gene[up$LS == x]
  if (length(g) < 5) next
  em <- tryCatch(enricher(g, TERM2GENE=BP, pvalueCutoff=1, qvalueCutoff=1, minGSSize=5, maxGSSize=500),
                 error=function(e) NULL)
  if (is.null(em) || !nrow(as.data.frame(em))) next
  d <- as.data.frame(em); d$LS <- x
  res[[x]] <- d[, c("LS","Description","p.adjust")]
}
R <- bind_rows(res)
R$term <- R$Description |> sub("^GOBP_","",x=_) |> gsub("_"," ",x=_) |> tolower()
# top terms per state, union
top <- R %>% group_by(LS) %>% slice_min(p.adjust, n=3, with_ties=FALSE) %>% ungroup()
terms <- unique(top$term)
if (length(terms) > 28) terms <- R %>% filter(term %in% terms) %>% group_by(term) %>%
  summarise(m=min(p.adjust)) %>% slice_min(m, n=28) %>% pull(term)
M <- R %>% filter(term %in% terms) %>% mutate(nl=-log10(p.adjust))
full <- expand.grid(term=terms, LS=ord$LS, stringsAsFactors=FALSE) %>%
  left_join(M[,c("LS","term","nl")], by=c("LS","term")) %>% mutate(nl=ifelse(is.na(nl),0,nl))
full <- left_join(full, ord[,c("LS","lab","group","dir")], by="LS")
# order columns fav-first (ord already is); order terms by the state where they peak
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
        plot.margin=margin(6,10,6,6))
wd <- max(6.5, 3.2 + 0.42*nrow(ord))
ggsave(file.path(OUT, paste0(OUTB,".pdf")), hp, width=wd, height=7.2)
ggsave(file.path(OUT, paste0(OUTB,".png")), hp, width=wd, height=7.2, dpi=200)
write.csv(R, file.path(SL, paste0("ls48/", OUTB, "_all.csv")), row.names=FALSE)
cat(sprintf("wrote %s: %d terms x %d states\n", OUTB, length(terms), nrow(ord)))
