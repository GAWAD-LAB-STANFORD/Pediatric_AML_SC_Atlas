#!/usr/bin/env Rscript
# LS10 deconvolved fraction (%) across TARGET cytogenetic subtypes, coloured by risk group.
# LS10 is zero-inflated (81% zero), so points are jittered and a diamond marks the mean;
# y is sqrt-scaled to show the spread. Shows LS10 enrichment in adverse-risk subtypes.
suppressMessages({library(ggplot2); library(scales)})
CODE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3"
source(file.path(CODE, "manuscript_palette.R"))
D3 <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
d <- read.csv(file.path(D3, "LS_survival.csv"), check.names=FALSE)
d <- d[!is.na(d$subtype2) & trimws(d$subtype2) != "", c("subtype2","LS_10")]; d$pct <- d$LS_10*100  # drop blank/unclassified cytogenetics
keep <- names(which(table(d$subtype2) >= 15)); d <- d[d$subtype2 %in% keep, ]
risk <- c("PML-RARA"="favorable","t(8;21)"="favorable","inv(16)"="favorable",
          "Normal"="intermediate","KMT2A(MLL)"="intermediate","Other"="intermediate",
          "mono7/del7q"="adverse","t(6;9)"="adverse","NUP98-r"="adverse","complex(>=3)"="adverse")
d$risk <- factor(risk[d$subtype2], levels=c("favorable","intermediate","adverse"))
ord <- names(sort(tapply(d$pct, d$subtype2, mean)))          # order by mean LS10 %
d$subtype2 <- factor(d$subtype2, levels=ord)
cols <- c(favorable=unname(PROG["favorable"]), intermediate=GREY, adverse=unname(PROG["poor"]))
g <- ggplot(d, aes(subtype2, pct, colour=risk)) +
  geom_jitter(width=0.18, height=0, size=0.7, alpha=0.4) +
  geom_boxplot(aes(fill=risk), outlier.shape=NA, width=0.55, alpha=0.25, linewidth=0.4, colour="grey40") +
  stat_summary(fun=mean, geom="point", shape=23, size=2.6, fill="white", colour="black") +
  scale_colour_manual(values=cols, name="risk", na.translate=FALSE) + scale_fill_manual(values=cols, guide="none", na.translate=FALSE) +
  scale_y_sqrt(breaks=c(0,0.5,1,2,4,8)) +
  labs(x="cytogenetic / molecular subtype", y="LS10 deconvolved fraction (%)") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=40, hjust=1), legend.position="top",
        plot.caption=element_text(size=8, hjust=0))
ggsave(file.path(OUT, "Figure_2_LS10_by_cyto.pdf"), g, width=6.4, height=4.2)
ggsave(file.path(OUT, "Figure_2_LS10_by_cyto.png"), g, width=6.4, height=4.2, dpi=200)
cat("wrote Figure_2_LS10_by_cyto\n")
print(round(sort(tapply(d$pct, d$subtype2, mean)),3))
