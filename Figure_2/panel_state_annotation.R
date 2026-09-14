#!/usr/bin/env Rscript
# S10: annotation of the 49 leukemic states (corrected scheme; Jaccard >= 0.35).
# Heatmap of curated hematopoietic lineage-program markers (z per gene across states).
# States labelled by DIFFERENTIATION PROGRAM ONLY (no cytogenetics). Right bars: cross-
# cytogenetic breadth (effective # subtypes) and patient-dominance class. Data =
# compute_state_annotation.py (current 49-state h5ad). Replaces the legacy 48-state figure.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
D2  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_2/v3"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/05_Supplementary_Figures"
PROG <- list("HSC/MPP"=c("CD34","HLF","AVP","CRHBP","HOXA9","MEIS1","MLLT3","PROM1","MECOM","GATA2"),
 "GMP"=c("MPO","ELANE","AZU1","PRTN3","CTSG","CEBPE","LYST","CFD"),
 "Monocyte"=c("LYZ","CD14","CSF1R","FCN1","VCAN","S100A8","S100A9","CD68"),
 "Macrophage"=c("C1QA","C1QB","C1QC","MRC1","SIGLEC1","TREM1"),
 "cDC"=c("CD1C","FCER1A","CLEC10A","CD207","CD1E","IRF8"),"pDC"=c("IRF7","LILRA4","GZMB","IL3RA"),
 "Megakaryocyte"=c("PPBP","PF4","ITGA2B","GP9","GP1BA","VWF","TUBB1","PF4V1"),
 "Erythroid"=c("GATA1","KLF1","HBB","HBA1","GYPA","ALAS2","TFRC","AHSP"),
 "Neutrophil"=c("FUT4","CEACAM8","MMP8","CD177","LTF","LCN2","CAMP","MMP9"),
 "B/plasma"=c("CD19","MS4A1","CD79A","CD79B","IGHM","JCHAIN","BLNK"),
 "T/NK"=c("CD3D","CD3E","IL7R","NKG7","GNLY","CD2","GZMH"),
 "proliferative"=c("MKI67","TOP2A","CDK1","UBE2C","CENPA","BIRC5","PCNA","NUSAP1","BUB1B"))
g2p <- setNames(rep(names(PROG),lengths(PROG)),unlist(PROG))
prog_ord <- c(names(PROG),"unresolved")
z   <- read.csv(file.path(D2,"state_annotation_heatmap_z.csv"),row.names=1,check.names=FALSE)
ann <- read.csv(file.path(D2,"state_annotation_meta.csv"),check.names=FALSE)
ann$lp <- factor(ann$lineage_program,levels=prog_ord)
ann <- ann[order(ann$lp,-ann$n_cells),]; state_ord <- ann$LS
ann$label <- paste0(sub("_","",ann$LS),"  ",ann$lineage_program)
lab <- setNames(ann$label,ann$LS)
mk <- colnames(z); mk <- mk[mk %in% names(g2p)]
long <- do.call(rbind,lapply(state_ord,function(s) data.frame(state=s,gene=mk,z=as.numeric(z[s,mk]),program=g2p[mk])))
long$state <- factor(long$state,levels=rev(state_ord),labels=rev(lab[state_ord]))
long$gene <- factor(long$gene,levels=mk); long$program <- factor(long$program,levels=names(PROG))
ph <- ggplot(long,aes(gene,state,fill=z))+geom_tile()+
  facet_grid(~program,scales="free_x",space="free_x",switch="x")+
  scale_fill_gradient2(low="#2166ac",mid="white",high="#b2182b",midpoint=0,limits=c(-2.5,2.5),oob=squish,name="z-score")+
  labs(x=NULL,y=NULL)+theme_minimal(base_size=7)+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=4.5),axis.text.y=element_text(size=6),
        panel.grid=element_blank(),strip.text.x=element_text(angle=90,size=6,face="bold"),
        strip.placement="outside",panel.spacing=unit(0.05,"lines"),legend.position="left")
ad <- ann; ad$state <- factor(ad$LS,levels=rev(state_ord),labels=rev(lab[state_ord]))
b1 <- ggplot(ad,aes("cyto\nbreadth",state,fill=eff_n_cyto))+geom_tile(colour="white",linewidth=0.3)+
  scale_fill_viridis_c(option="G",name="eff. #\nsubtypes")+
  labs(x=NULL,y=NULL)+theme_minimal(base_size=7)+
  theme(axis.text.y=element_blank(),axis.text.x=element_text(size=6),panel.grid=element_blank(),legend.position="right")
pcol <- c(`multi-pt`="#08519c",`few-pt`="#6baed6",`single-pt`="#c6dbef")
ad$patient_class <- factor(ad$patient_class,levels=names(pcol))
b2 <- ggplot(ad,aes("patients",state,fill=patient_class))+geom_tile(colour="white",linewidth=0.3)+
  scale_fill_manual(values=pcol,name=NULL,drop=FALSE)+labs(x=NULL,y=NULL)+theme_minimal(base_size=7)+
  theme(axis.text.y=element_blank(),axis.text.x=element_text(size=6),panel.grid=element_blank(),legend.position="right")
fig <- ph + b1 + b2 + plot_layout(widths=c(1,0.04,0.04))
ggsave(file.path(OUT,"Figure_S10__leukemic_state_annotation.pdf"),fig,width=15,height=10)
ggsave(file.path(OUT,"Figure_S10__leukemic_state_annotation.png"),fig,width=15,height=10,dpi=200)
cat("wrote 49-state annotation heatmap (S10)\n")
