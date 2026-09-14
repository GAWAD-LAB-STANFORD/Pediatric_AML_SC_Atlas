#!/usr/bin/env Rscript
# Figure 3 panel E: GO:BP over-representation (clusterProfiler enricher vs MSigDB
# C5:BP) of the up-regulated genes per CBF-prognostic state; heatmap of -log10 q.
suppressPackageStartupMessages({library(clusterProfiler); library(msigdbr); library(dplyr); library(ggplot2); library(scales); library(stringr)})
source(file.path(dirname(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))),"config.R"))
get_bp <- function() for (call in list(
    function() msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP"),
    function() msigdbr(species="Homo sapiens", category="C5", subcategory="GO:BP")))
  { r<-tryCatch(call(),error=function(e)NULL); if(!is.null(r)&&nrow(r)) return(r) }
BP <- get_bp() %>% dplyr::select(gs_name, gene_symbol)
up  <- read.csv(F3$UPGENES); ord <- read.csv(file.path(F3$DATA,"LS_cbf_order.csv"))
res <- list()
for (i in seq_len(nrow(ord))) { x<-ord$LS[i]; g<-up$gene[up$LS==x]; if(length(g)<5) next
  em<-tryCatch(enricher(g,TERM2GENE=BP,pvalueCutoff=1,qvalueCutoff=1,minGSSize=5,maxGSSize=500),error=function(e)NULL)
  if(is.null(em)||!nrow(as.data.frame(em))) next
  d<-as.data.frame(em); d$LS<-x; res[[x]]<-d[,c("LS","Description","p.adjust")] }
R<-bind_rows(res); R$term<-R$Description|>sub("^GOBP_","",x=_)|>gsub("_"," ",x=_)|>tolower()
top<-R%>%group_by(LS)%>%slice_min(p.adjust,n=3,with_ties=FALSE)%>%ungroup(); terms<-unique(top$term)
if(length(terms)>28) terms<-R%>%filter(term%in%terms)%>%group_by(term)%>%summarise(m=min(p.adjust))%>%slice_min(m,n=28)%>%pull(term)
M<-R%>%filter(term%in%terms)%>%mutate(nl=-log10(p.adjust))
full<-expand.grid(term=terms,LS=ord$LS,stringsAsFactors=FALSE)%>%left_join(M[,c("LS","term","nl")],by=c("LS","term"))%>%
  mutate(nl=ifelse(is.na(nl),0,nl))%>%left_join(ord[,c("LS","lab","dir")],by="LS")
peak<-full%>%group_by(term)%>%slice_max(nl,n=1,with_ties=FALSE)%>%mutate(ci=match(LS,ord$LS))%>%select(term,ci)%>%arrange(ci)
full$term<-factor(full$term,levels=rev(peak$term)); full$lab<-factor(full$lab,levels=ord$lab)
full$term<-factor(full$term,levels=levels(full$term),labels=str_wrap(levels(full$term),34))
labcol<-ifelse(ord$dir=="favorable",unname(PROG["favorable"]),unname(PROG["poor"]))
hp<-ggplot(full,aes(lab,term,fill=nl))+geom_tile(colour="white",linewidth=0.4)+
  geom_vline(xintercept=sum(ord$dir=="favorable")+0.5,colour="grey30",linewidth=0.6)+
  scale_fill_seq(limits=c(0,6),oob=squish,name=expression(-log[10]~q))+labs(x=NULL,y=NULL)+
  theme_minimal(base_size=9)+theme(axis.text.x=element_text(angle=45,hjust=1,face="bold",size=8.5,colour=labcol),
    axis.text.y=element_text(size=7.5),panel.grid=element_blank(),plot.margin=margin(6,10,6,6))
wd<-max(6.5,3.2+0.42*nrow(ord))
ggsave(file.path(F3$OUT,"panelE_cbf_go.pdf"),hp,width=wd,height=7.2)
ggsave(file.path(F3$OUT,"panelE_cbf_go.png"),hp,width=wd,height=7.2,dpi=200)
cat(sprintf("wrote panelE_cbf_go: %d terms x %d states\n",length(terms),nrow(ord)))
