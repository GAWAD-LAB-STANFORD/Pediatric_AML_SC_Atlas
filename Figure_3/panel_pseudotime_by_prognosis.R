#!/usr/bin/env Rscript
# Pseudotime ordering (dpt AND Monocle3) with leukemic states coloured by whole-cohort EFS prognosis
# (favorable / poor / n.s.) and normal compartments as bold purple reference anchors. Data =
# pseudotime_by_group.csv, monocle_pseudotime.csv, LS_wholecohort_prognosis.csv. No fabrication.
suppressMessages({library(ggplot2); library(dplyr)})
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
D2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_3_panels"
pr  <- read.csv(file.path(D2,"LS_wholecohort_prognosis.csv"))[,c("LS","prognosis")]
pcol <- c(favorable="#00798c", poor="#d1495b", `n.s.`="#9aa0a6", normal="#7b3294")
mkfig <- function(csv, val, xlab, outfile){
  d <- read.csv(file.path(D3,csv)); d <- d[is.finite(d[[val]]),]
  q <- d %>% group_by(group) %>% summarise(med=median(.data[[val]]), lo=quantile(.data[[val]],.25),
                                            hi=quantile(.data[[val]],.75), .groups="drop")
  q$normal <- grepl("^NORM_", q$group)
  q <- left_join(q, pr, by=c("group"="LS"))
  q$prognosis[q$normal] <- "normal"
  q$prognosis[is.na(q$prognosis)] <- "n.s."
  q$prognosis <- factor(ifelse(q$prognosis=="ns","n.s.",q$prognosis), levels=names(pcol))
  q$label <- ifelse(q$normal, paste0("norm ", sub("NORM_","",q$group)), sub("_","",q$group))
  q <- q %>% arrange(med); q$label <- factor(q$label, levels=q$label)
  p <- ggplot(q, aes(med, label, colour=prognosis)) +
    geom_segment(aes(x=lo,xend=hi,y=label,yend=label), linewidth=0.5) +
    geom_point(aes(shape=normal, size=normal)) +
    scale_colour_manual(values=pcol, name="EFS prognosis") +
    scale_shape_manual(values=c(`FALSE`=16,`TRUE`=18), guide="none") +
    scale_size_manual(values=c(`FALSE`=2,`TRUE`=3.4), guide="none") +
    labs(x=xlab, y=NULL) + theme_minimal(base_size=9) +
    theme(panel.grid.major.y=element_line(colour="grey93"), panel.grid.minor=element_blank(),
          axis.text.y=element_text(size=6.4), legend.position="right", legend.key.size=unit(0.32,"cm"))
  ggsave(file.path(OUT,paste0(outfile,".pdf")), p, width=7, height=11)
  ggsave(file.path(OUT,paste0(outfile,".png")), p, width=7, height=11, dpi=170)
}
mkfig("pseudotime_by_group.csv","dpt","diffusion pseudotime (dpt, rooted at normal HSPC)","Figure_3_dpt_ordering_prognosis")
mkfig("monocle_pseudotime.csv","monocle_pt","Monocle3 pseudotime (rooted at normal HSPC)","Figure_3_monocle_ordering_prognosis")
cat("wrote dpt + Monocle3 ordering figures coloured by prognosis\n")
