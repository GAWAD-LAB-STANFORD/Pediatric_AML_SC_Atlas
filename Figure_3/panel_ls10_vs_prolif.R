#!/usr/bin/env Rscript
# Figure 3C (rebuilt): what makes LS10 unique AMONG proliferating blasts. Differential expression
# of LS10 vs the other three proliferation-high states (LS3, LS20, LS39), which cancels out
# proliferation itself. Top 10 up and top 10 down genes by log2 fold-change (Wilcoxon, BH-adjusted
# p<0.05, expressed in >25% of the higher group). UP = an E2F-driven S-phase/replication program;
# DOWN = a myeloid-differentiation program, i.e. LS10 couples high replication to a maturation
# block. Data = compute (LS10_vs_prolif_DE.csv). No fabricated values.
suppressMessages({library(ggplot2); library(dplyr)})
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_2_panels"
de <- read.csv(file.path(D3,"LS10_vs_prolif_DE.csv"))
up <- de %>% filter(pvals_adj<0.05, logfoldchanges> 0.5, pct_LS10 >0.25) %>% arrange(desc(logfoldchanges)) %>% head(10)
dn <- de %>% filter(pvals_adj<0.05, logfoldchanges< -0.5, pct_other>0.25) %>% arrange(logfoldchanges) %>% head(10)
d <- bind_rows(up %>% mutate(dir="up in LS10"), dn %>% mutate(dir="down in LS10"))
d$names <- factor(d$names, levels=d$names[order(d$logfoldchanges)])
d$dir <- factor(d$dir, levels=c("up in LS10","down in LS10"))
p <- ggplot(d, aes(logfoldchanges, names, fill=dir)) +
  geom_col(width=0.72) + geom_vline(xintercept=0, colour="grey50", linewidth=0.4) +
  scale_fill_manual(values=c("up in LS10"="#d1495b","down in LS10"="#00798c"), name=NULL) +
  labs(x="log2 fold-change (LS10 vs LS3/LS20/LS39)", y=NULL,
       title="What distinguishes LS10 among proliferating states",
       subtitle="up: E2F/S-phase replication   |   down: myeloid differentiation") +
  theme_minimal(base_size=11) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=8.5), plot.title=element_text(size=11,face="bold"),
        plot.subtitle=element_text(size=8.5,colour="grey30"), legend.position=c(0.82,0.15),
        legend.background=element_rect(fill="white",colour=NA))
ggsave(file.path(OUT,"Figure_3_LS10_vs_prolif.pdf"), p, width=5.6, height=4.6)
ggsave(file.path(OUT,"Figure_3_LS10_vs_prolif.png"), p, width=5.6, height=4.6, dpi=200)
cat("wrote 3C LS10-vs-proliferation DE: up=",paste(up$names,collapse=","),"| down=",paste(dn$names,collapse=","),"\n")
