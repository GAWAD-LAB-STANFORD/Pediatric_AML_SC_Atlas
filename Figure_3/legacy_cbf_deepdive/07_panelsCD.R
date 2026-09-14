#!/usr/bin/env Rscript
# Figure 3 panels C (focus-state: distinct + targetable) and D (HOX expression).
suppressMessages({library(ggplot2); library(scales)})
source(file.path(dirname(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))),"config.R"))
focus <- readLines(file.path(F3$DATA,"LS_focus_state.txt"))[1]; flab <- gsub("_","",focus)

## ---- Panel C: focus-state identity + targetable surface antigens ----
d <- read.csv(file.path(F3$DATA,"panelC_targetable.csv"))
cols <- c("Normal_HSPC","Favorable_mean","Poor_mean","Focus")
labs <- c("Normal\nHSPC","Favorable\n(mean)","Poor\n(mean)",flab)
long <- do.call(rbind, lapply(seq_along(cols), function(i) data.frame(gene=d$gene,group=d$group,col=labs[i],v=d[[cols[i]]])))
grporder <- c("not LT-HSC","focus signature","shared progenitor","not GMP (granule)","cycling","surface targets")
long$group<-factor(long$group,levels=grporder); long$gene<-factor(long$gene,levels=rev(d$gene)); long$col<-factor(long$col,levels=labs)
pc <- ggplot(long,aes(col,gene,fill=v))+geom_tile(colour="white",linewidth=0.4)+
  facet_grid(group~.,scales="free_y",space="free_y",switch="y")+
  scale_fill_seq(name="mean expr\nlog1p(CPM)")+labs(x=NULL,y=NULL)+theme_minimal(base_size=10)+
  theme(axis.text.x=element_text(size=8.5,face=ifelse(labs==flab,"bold","plain")),axis.text.y=element_text(size=8),
    panel.grid=element_blank(),strip.text.y.left=element_text(angle=0,hjust=1,size=7.5,face="bold"),
    strip.placement="outside",panel.spacing=unit(2,"pt"),plot.margin=margin(6,8,6,6))
ggsave(file.path(F3$OUT,"panelC_targetable.pdf"),pc,width=5.2,height=7.6)
ggsave(file.path(F3$OUT,"panelC_targetable.png"),pc,width=5.2,height=7.6,dpi=200)

## ---- Panel D: HOX expression across CBF-prognostic states ----
E <- read.csv(file.path(F3$DATA,"panelD_hox.csv"),row.names=1,check.names=FALSE); genes<-rownames(E)
ord <- read.csv(file.path(F3$DATA,"LS_cbf_order.csv")); ann<-read.csv(F3$ANNOT); prog<-setNames(ann$lineage_program,ann$LS)
stlev<-ord$LS; nfav<-sum(ord$dir=="favorable")
lab<-ifelse(prog[stlev]%in%c("unresolved",NA),gsub("_","",stlev),paste0(gsub("_","",stlev)," ",prog[stlev]))
labcol<-ifelse(ord$dir=="favorable",unname(PROG["favorable"]),unname(PROG["poor"]))
long2<-do.call(rbind,lapply(genes,function(g) data.frame(gene=g,LS=stlev,v=as.numeric(E[g,stlev]))))
long2$LS<-factor(long2$LS,levels=stlev,labels=lab); long2$gene<-factor(long2$gene,levels=rev(genes))
pd<-ggplot(long2,aes(LS,gene,fill=v))+geom_tile(colour="white",linewidth=0.3)+
  geom_vline(xintercept=nfav+0.5,colour="grey30",linewidth=0.6)+scale_fill_seq(name="mean expr\nlog1p(CPM)")+
  labs(x=NULL,y=NULL)+theme_minimal(base_size=10)+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8.5,face="bold",colour=labcol),
    axis.text.y=element_text(size=7.5),panel.grid=element_blank(),plot.margin=margin(6,10,6,6))
ggsave(file.path(F3$OUT,"panelD_hox.pdf"),pd,width=8.4,height=7.2)
ggsave(file.path(F3$OUT,"panelD_hox.png"),pd,width=8.4,height=7.2,dpi=200)
cat("wrote panelC_targetable and panelD_hox (focus =",focus,")\n")
