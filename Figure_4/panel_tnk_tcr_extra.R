#!/usr/bin/env Rscript
# Figure 4 extended TCR panels. Supplement Figure S12 (composite_supp_tcr.py) uses THREE of these:
#   A = (2) clonotype sharing heatmap (Jaccard; TCR as a lineage barcode)
#   B = (3) exhaustion score, expanded vs unexpanded clones (honest null)
#   C = (4) clonality by cytogenetic subtype (honest null; small per-sample n shown as n= labels)
# Interpretive text lives in the S12 legend, not on the panels (no on-figure subtitles). PDF + PNG saved.
# (1) per-subset downsampled clonality and (5) TRBV usage AML-vs-HBM are also built here but are NOT in
# the deployed supplement: (1) duplicates main Fig 4F (per-patient clonality); (5) was dropped because the
# HBM normal reference is only 2 donors (underpowered comparator).
suppressMessages({library(ggplot2); library(dplyr); library(scales)})
source("/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/07_Code/Figure_3/manuscript_palette.R")
V2 <-"/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_4/v2"
OUT<-"/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_4_panels"
sv <- read.csv(file.path(V2,"tnk_subset_survival.csv")) %>% filter(outcome=="EFS") %>%
      mutate(dir=ifelse(FDR<0.10, ifelse(HR>1,"poor","favorable"),"n.s.")) %>% select(subset,dir)
PC <- c(favorable=unname(PROG["favorable"]),poor=unname(PROG["poor"]),`n.s.`="#9aa0a6")

## (1) Shannon clonality (downsampled) per subset
d1 <- merge(read.csv(file.path(V2,"tcr_clonality_downsampled.csv")), sv, by="subset", all.x=TRUE)
d1$dir[is.na(d1$dir)]<-"n.s."; d1 <- d1 %>% arrange(clonality_ds); d1$subset<-factor(d1$subset,levels=d1$subset)
g1 <- ggplot(d1, aes(clonality_ds, subset, fill=factor(dir,levels=names(PC)))) +
  geom_col(width=.7,colour="grey30",linewidth=.2) +
  geom_errorbarh(aes(xmin=pmax(0,clonality_ds-clonality_sd), xmax=clonality_ds+clonality_sd), height=.25, linewidth=.4) +
  scale_fill_manual(values=PC,name="EFS prognosis") + scale_x_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(x="clonality (1 - normalized Shannon entropy), downsampled to n=100", y=NULL) +
  theme_classic(base_size=11)+theme(legend.position=c(.98,.05),legend.justification=c(1,0))
ggsave(file.path(OUT,"Figure_4_tcr_clonality.png"),g1,width=5.8,height=4.2,dpi=200)

## (2) clonotype sharing heatmap (Jaccard)
sh <- read.csv(file.path(V2,"tcr_sharing.csv"))
ORD <- c("Naïve CD4 T","Memory CD4 T","Treg","Naïve CD8 T","Memory CD8 T","GZMK CD8 T","GZMB CD8 T","MAIT")
sh <- sh %>% filter(subset_a%in%ORD, subset_b%in%ORD)
sh$diag <- as.character(sh$subset_a)==as.character(sh$subset_b)
sh$jaccard[sh$diag] <- NA                      # blank the self-overlap diagonal so the scale reflects between-subset sharing
sh$subset_a<-factor(sh$subset_a,levels=ORD); sh$subset_b<-factor(sh$subset_b,levels=rev(ORD))
sh$lab <- ifelse(sh$n_shared>0 & !sh$diag, sh$n_shared, "")
g2 <- ggplot(sh, aes(subset_a, subset_b, fill=jaccard)) + geom_tile(colour="grey85") +
  geom_text(aes(label=lab), size=2.8, colour="grey15") +
  scale_fill_seq(name="Jaccard\n(clonotype overlap)", na.value="grey92") +
  labs(x=NULL,y=NULL) +
  theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())
ggsave(file.path(OUT,"Figure_4_tcr_sharing.png"),g2,width=5.6,height=4.8,dpi=200)
ggsave(file.path(OUT,"Figure_4_tcr_sharing.pdf"),g2,width=5.6,height=4.8)

## (3) exhaustion score: expanded vs unexpanded (cytotoxic/memory CD8)
ex <- read.csv(file.path(V2,"tcr_exhaustion_by_expansion.csv"))
ex <- ex %>% filter(subset %in% c("GZMK CD8 T","GZMB CD8 T","Memory CD8 T")) %>%
      tidyr::pivot_longer(c(expanded,unexpanded),names_to="clone",values_to="exh")
g3 <- ggplot(ex, aes(subset, exh, fill=clone)) +
  geom_col(position="dodge", width=.7, colour="grey30", linewidth=.2) +
  scale_fill_manual(values=c(expanded="#AE2012",unexpanded="#0A9396"),name=NULL) +
  labs(x=NULL,y="exhaustion score (PDCD1/HAVCR2/LAG3/TOX)") +
  theme_classic(base_size=11)+theme(axis.text.x=element_text(angle=20,hjust=1))
ggsave(file.path(OUT,"Figure_4_tcr_exhaustion.png"),g3,width=4.8,height=4.2,dpi=200)
ggsave(file.path(OUT,"Figure_4_tcr_exhaustion.pdf"),g3,width=4.8,height=4.2)

## (4) clonality by cytogenetic subtype (per AML sample)
cy <- read.csv(file.path(V2,"tcr_clonality_by_cyto.csv"))
cy$Cytogenetic <- reorder(cy$Cytogenetic, cy$clonality, median)
nlab <- cy %>% dplyr::count(Cytogenetic)                       # per-subtype sample n (small; show honestly)
g4 <- ggplot(cy, aes(Cytogenetic, clonality)) +
  stat_summary(fun=median, fun.min=median, fun.max=median, geom="crossbar", width=.5, colour="grey45", linewidth=.35) +
  geom_jitter(width=.12,size=1.8,colour="#AE2012",alpha=.8) +
  geom_text(data=nlab, aes(Cytogenetic, y=-Inf, label=paste0("n=",n)), vjust=-0.6, size=2.4, colour="grey40", inherit.aes=FALSE) +
  labs(x="cytogenetic / molecular subtype", y="T-cell clonality (per sample); grey bar = median") +
  theme_classic(base_size=10)+theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(OUT,"Figure_4_tcr_clonality_cyto.png"),g4,width=6.2,height=4.2,dpi=200)
ggsave(file.path(OUT,"Figure_4_tcr_clonality_cyto.pdf"),g4,width=6.2,height=4.2)

## (5) TRBV usage AML vs HBM (top V genes)
vb <- read.csv(file.path(V2,"tcr_trbv_usage.csv"))
top <- vb %>% group_by(Vgene) %>% summarise(m=sum(frac)) %>% slice_max(m,n=18) %>% pull(Vgene)
vb <- vb %>% filter(Vgene%in%top); vb$Vgene <- reorder(vb$Vgene, vb$frac, max)
g5 <- ggplot(vb, aes(frac*100, Vgene, fill=group)) +
  geom_col(position="dodge", width=.7, colour="grey30", linewidth=.2) +
  scale_fill_manual(values=c(AML=unname(PROG["poor"]),HBM=unname(PROG["favorable"])),name=NULL) +
  labs(x="% of TCR+ cells", y=NULL, title="TRBV usage") +
  theme_classic(base_size=10)+theme(legend.position=c(.98,.05),legend.justification=c(1,0),axis.text.y=element_text(size=7))
ggsave(file.path(OUT,"Figure_4_tcr_trbv.png"),g5,width=5.2,height=5.2,dpi=200)
cat("wrote 5 TCR panels: clonality, sharing, exhaustion, clonality_cyto, trbv\n")
