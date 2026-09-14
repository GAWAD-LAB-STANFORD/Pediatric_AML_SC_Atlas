#!/usr/bin/env Rscript
# Figure 2C (v3, corrected 49 states): 14 whole-cohort EFS-prognostic states x
# {cytogenetic | molecular} subtypes (dot size = mean deconvolved %, fill = -log10 FDR
# of enrichment) + a stemness strip. Enrichment is computed here from the corrected
# CIBERSORTx fractions (LS_survival.csv) + per-patient molecular status. Single-patient
# states are marked with *. Harmonized ltc palette. No fabricated values.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_2"
source(file.path(CODE, "manuscript_palette.R"))
D2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
MOL <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE/08_Source_Data/Figure_2/Figure2CD_TARGET_PAC_perpatient.csv"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# ---- data ----
sv  <- read.csv(file.path(D3, "LS_survival.csv"), check.names=FALSE)
lsN <- grep("^LS_", names(sv), value=TRUE)
mol <- read.csv(MOL, check.names=FALSE)
recode <- function(x){ x <- ifelse(x %in% c("NO"), "No", x); ifelse(x %in% c("Yes","No"), x, NA) }
molmap <- data.frame(sample=mol$sample,
  `FLT3-ITD`=recode(mol[["FLT3/ITD positive?"]]), `NPM1-mut`=recode(mol[["NPM mutation"]]),
  `CEBPA-mut`=recode(mol[["CEBPA mutation"]]), `WT1-mut`=recode(mol[["WT1 mutation"]]), check.names=FALSE)
sv <- merge(sv, molmap, by="sample", all.x=TRUE)

# ---- enrichment: mean % + one-sided Wilcoxon (in-group vs rest), BH-FDR across all tests ----
enrich <- function(dat, groupcol, groups){
  r <- data.frame()
  for (g in groups){
    inx <- which(dat[[groupcol]] == g); if (length(inx) < 15) next
    for (ls in lsN){ fr <- dat[[ls]]
      p <- tryCatch(wilcox.test(fr[inx], fr[-inx], alternative="greater")$p.value, error=function(e) NA)
      r <- rbind(r, data.frame(subtype=g, n=length(inx), LS=ls, mean_pct=mean(fr[inx])*100, pval=p)) }
  }
  r$fdr <- p.adjust(r$pval, "BH"); r$neglog10fdr <- -log10(r$fdr); r
}
cyto_lv <- c("PML-RARA","t(8;21)","inv(16)","Normal","KMT2A(MLL)","NUP98-r","t(6;9)","mono7/del7q","complex(>=3)","Other")
mol_lv  <- c("FLT3-ITD","NPM1-mut","CEBPA-mut","WT1-mut")
dc_all <- enrich(sv, "subtype2", cyto_lv)
dm_all <- do.call(rbind, lapply(mol_lv, function(mc){
  d2 <- sv[!is.na(sv[[mc]]), ]; d2$grp <- d2[[mc]]
  inx <- which(d2$grp=="Yes"); if (length(inx) < 15) return(NULL)
  do.call(rbind, lapply(lsN, function(ls){ fr <- d2[[ls]]
    data.frame(subtype=mc, n=length(inx), LS=ls, mean_pct=mean(fr[inx])*100,
      pval=tryCatch(wilcox.test(fr[inx], fr[-inx], alternative="greater")$p.value, error=function(e) NA)) })) }))
dm_all$fdr <- p.adjust(dm_all$pval, "BH"); dm_all$neglog10fdr <- -log10(dm_all$fdr)
write.csv(dc_all, file.path(D2, "LS_cyto_enrichment.csv"), row.names=FALSE)
write.csv(dm_all, file.path(D2, "LS_mol_enrichment.csv"),  row.names=FALSE)

# ---- prognostic states + labels ----
pr <- read.csv(file.path(D2, "LS_wholecohort_prognosis.csv")); pr <- pr[pr$prognosis %in% c("favorable","poor"), ]
pr$LSc  <- gsub("_","",pr$LS)
pr$star <- ifelse(pr$n_patients <= 1, "*", "")
pr$progn <- ifelse(pr$lineage_program=="unresolved", "", pr$lineage_program)
pr$lab  <- ifelse(pr$progn=="", paste0(pr$LSc,pr$star), paste0(pr$LSc,pr$star,"  ",pr$progn))
pr <- pr[order(pr$prognosis=="poor", pr$p), ]; pr$lab <- make.unique(pr$lab)
lv <- pr$lab; nfav <- sum(pr$prognosis=="favorable"); yline <- nrow(pr) - nfav + 0.5

mkdat <- function(d, xlv){ d <- d[d$LS %in% pr$LS, ]; d <- merge(d, pr[,c("LS","lab")], by="LS")
  d$lab <- factor(d$lab, levels=rev(lv)); d$subtype <- factor(d$subtype, levels=xlv); d }
dc <- mkdat(dc_all, cyto_lv); dm <- mkdat(dm_all, mol_lv)
szmax <- ceiling(max(c(dc$mean_pct, dm$mean_pct), na.rm=TRUE))

dotlayer <- function(d, ylabs, ttl){
  ggplot(d, aes(subtype, lab, size=mean_pct, fill=neglog10fdr)) +
    geom_hline(yintercept=yline, colour="grey40", linewidth=0.5) +
    geom_point(shape=21, colour="grey30", stroke=0.3) +
    scale_fill_seq(limits=c(0,10), oob=squish, name=expression(-log[10]~FDR)) +
    scale_size(range=c(0.4,8.5), limits=c(0,szmax), name="mean % (deconv.)") +
    labs(x=NULL, y=NULL) + ggtitle(ttl) + theme_minimal(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1, face="bold"),
          axis.text.y=if (ylabs) element_text(size=9) else element_blank(),
          plot.title=element_text(size=10, face="plain"),
          panel.grid.major=element_line(colour="grey93")) }

blk <- data.frame(y=c(nrow(pr)-nfav/2+0.5, (nrow(pr)-nfav)/2+0.5), lab=c("Favorable","Poor"))
pblk <- ggplot(blk, aes(1, y, label=lab, colour=lab)) +
  geom_text(angle=90, fontface="bold", size=4.4) +
  scale_colour_manual(values=c(Favorable=unname(PROG["favorable"]), Poor=unname(PROG["poor"])), guide="none") +
  scale_y_continuous(limits=c(0.5, nrow(pr)+0.5), expand=c(0,0)) +
  scale_x_continuous(limits=c(0.5,1.5)) + theme_void()
pc <- dotlayer(dc, TRUE, "Cytogenetic"); pm <- dotlayer(dm, FALSE, "Molecular")
prs <- pr; prs$lab <- factor(prs$lab, levels=rev(lv))
ps <- ggplot(prs, aes("HSC/LSC", lab, fill=stem_z)) +
  geom_hline(yintercept=yline, colour="grey40", linewidth=0.5) +
  geom_tile(colour="white", linewidth=0.4) +
  scale_fill_div(midpoint=0, limits=c(-1,1), oob=squish, name="stemness\n(z)") +
  labs(x=NULL, y=NULL) + theme_minimal(base_size=10) +
  theme(axis.text.y=element_blank(), axis.text.x=element_text(angle=45, hjust=1, size=9), panel.grid=element_blank())
fig <- pblk + pc + pm + ps + plot_layout(widths=c(0.05,1,0.5,0.10), guides="collect") +
  plot_annotation(theme=theme(plot.margin=margin(12,8,6,8)))  # "* = single-patient state" defined in the figure legend
ggsave(file.path(OUT, "Figure_2C_states_subtypes.pdf"), fig, width=12, height=5.6)
ggsave(file.path(OUT, "Figure_2C_states_subtypes.png"), fig, width=12, height=5.6, dpi=200)
cat(sprintf("wrote Figure_2C: %d prognostic states, %d cyto x %d mol tests; size max %.0f%%\n",
    nrow(pr), nrow(dc_all), nrow(dm_all), szmax))
