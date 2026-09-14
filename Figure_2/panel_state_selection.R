#!/usr/bin/env Rscript
# S3: selection of the 49 leukemic states (corrected scheme). Number of reproducible states vs
# the per-cluster Jaccard reproducibility cutoff (Leiden res 4.5). The chosen cutoff 0.35 retains
# 49 reproducible states, keeping marginal but reproducible cross-cytogenetic states (incl. the
# smallest, e.g. the APL/PML-RARA state); a strict 0.60 would keep 22. Data = current v3
# jaccard_percluster.csv (config.stable_raws uses the same >=0.35 rule). Replaces the legacy
# 48-state figure, which read a superseded jaccard table. No fabricated values.
suppressMessages({library(ggplot2)})
JPC <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3/jaccard_percluster.csv"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/05_Supplementary_Figures"
jac <- read.csv(JPC); jac <- jac[jac$res==4.5,]
thr <- seq(0.30,0.75,0.05)
d <- data.frame(thr=thr, n=sapply(thr, function(t) sum(jac$mean_jaccard>=t)))
nchosen <- d$n[d$thr==0.35]; nstrict <- d$n[d$thr==0.60]
p <- ggplot(d,aes(thr,n))+
  geom_vline(xintercept=0.35,linetype="dashed",colour="#d1495b")+
  geom_vline(xintercept=0.60,linetype="dotted",colour="grey60")+
  geom_line(colour="#00798c",linewidth=0.9)+geom_point(colour="#00798c",size=2.6)+
  annotate("text",x=0.355,y=nchosen+3,label=paste0("chosen 0.35\n(",nchosen," states)"),colour="#d1495b",size=3.4,hjust=0)+
  annotate("text",x=0.60,y=nstrict-4,label=paste0("strict 0.60\n(",nstrict,")"),colour="grey45",size=3.2)+
  scale_x_continuous(breaks=thr)+
  labs(x="per-cluster Jaccard reproducibility cutoff (Leiden res 4.5)",
       y="number of reproducible leukemic states",
       title="Selection of the 49 reproducible leukemic states")+
  theme_minimal(base_size=12)+
  theme(panel.grid.minor=element_blank(),plot.title=element_text(size=12,face="bold"))
ggsave(file.path(OUT,"Figure_S3__leukemic_state_selection.pdf"),p,width=6.4,height=4.4)
ggsave(file.path(OUT,"Figure_S3__leukemic_state_selection.png"),p,width=6.4,height=4.4,dpi=200)
cat(sprintf("wrote 49-state selection figure (S3): chosen=%d, strict=%d\n",nchosen,nstrict))
