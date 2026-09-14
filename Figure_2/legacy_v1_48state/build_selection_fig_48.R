#!/usr/bin/env Rscript
# Supplemental: selection of the 48 leukemic states.
# A: number of reproducible states vs the Jaccard cutoff (res 4.5); 0.35 chosen
#    (48 states) to retain marginal cross-cytogenetic states incl. the APL state.
# B: bulk-deconvolution separability vs cluster number (why bulk resolves a subset).
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
SC  <- "/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/05_Supplementary_Figures/Figure_S_states"
jac <- read.csv(file.path(SC,"jaccard_hi/jaccard_percluster.csv")); jac <- jac[jac$res==4.5,]
thr <- seq(0.30,0.75,0.05)
nst <- sapply(thr,function(t) sum(jac$mean_jaccard>=t))
d <- data.frame(thr=thr,n=nst)
pA <- ggplot(d,aes(thr,n))+
  geom_vline(xintercept=0.35,linetype="dashed",colour="#d1495b")+
  geom_vline(xintercept=0.60,linetype="dotted",colour="grey60")+
  geom_line(colour="#00798c",linewidth=0.8)+geom_point(colour="#00798c",size=2.2)+
  annotate("text",x=0.35,y=52,label="chosen 0.35\n(48 states)",colour="#d1495b",size=3,hjust=0.1)+
  annotate("text",x=0.60,y=26,label="strict 0.6\n(30)",colour="grey50",size=3)+
  scale_x_continuous(breaks=thr)+
  labs(x="per-cluster Jaccard reproducibility cutoff",y="number of reproducible states",
       title="A  State count vs reproducibility cutoff (res 4.5)")+
  theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),plot.title=element_text(size=11,face="bold"))
sep <- read.csv(file.path(SC,"signatures_all/separability_all.csv"))
sc_ <- max(sep$cond_num)/max(sep$`pairs_r..95`+1)
pB <- ggplot(sep,aes(K))+
  geom_vline(xintercept=48,linetype="dashed",colour="#d1495b")+
  geom_line(aes(y=cond_num),colour="#00798c",linewidth=0.8)+geom_point(aes(y=cond_num),colour="#00798c",size=1.8)+
  geom_line(aes(y=`pairs_r..95`*sc_),colour="#d1495b",linewidth=0.8,linetype="dashed")+
  geom_point(aes(y=`pairs_r..95`*sc_),colour="#d1495b",size=1.8,shape=15)+
  scale_y_continuous(name="signature condition number",sec.axis=sec_axis(~ ./sc_,name="collinear signature pairs (r > 0.95)"))+
  labs(x="number of clusters",title="B  Bulk-deconvolution separability")+
  theme_minimal(base_size=11)+theme(axis.title.y=element_text(colour="#00798c"),axis.title.y.right=element_text(colour="#d1495b"),
    panel.grid.minor=element_blank(),plot.title=element_text(size=11,face="bold"))
ggsave(file.path(OUT,"Figure_S_state_selection.pdf"),pA+pB+plot_layout(widths=c(1.1,1)),width=11,height=4.3)
ggsave(file.path(OUT,"Figure_S_state_selection.png"),pA+pB+plot_layout(widths=c(1.1,1)),width=11,height=4.3,dpi=200)
cat("wrote 48-state selection figure\n")
