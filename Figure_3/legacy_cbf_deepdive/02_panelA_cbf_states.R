#!/usr/bin/env Rscript
# Figure 3 panel A: CBF-prognostic states x CBF cytogenetic enrichment + stemness
# strip + within-CBF EFS Cox HR(p). States/HR/p come from 01_cbf_prognosis.R.
suppressMessages({library(ggplot2); library(patchwork); library(scales)})
source(file.path(dirname(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE))),"config.R"))
ord <- read.csv(file.path(F3$DATA,"LS_cbf_order.csv"))
d   <- read.csv(F3$SURV, check.names=FALSE); cbf <- d[d$subtype2 %in% F3$CBF, ]
ann <- read.csv(F3$ANNOT); stemv <- setNames(ann$stem_z, ann$LS)
lv <- ord$lab; nfav <- sum(ord$dir=="favorable"); yline <- nrow(ord)-nfav+0.5
labcol <- ifelse(ord$dir=="favorable", unname(PROG["favorable"]), unname(PROG["poor"]))
# one-sided MWU enrichment of each state in each CBF subtype
enr <- do.call(rbind, lapply(seq_len(nrow(ord)), function(i){ x<-ord$LS[i]; fr<-cbf[[x]]; g<-cbf$subtype2
  do.call(rbind, lapply(F3$CBF, function(L){ a<-fr[g==L]; b<-fr[g!=L]
    p<-tryCatch(wilcox.test(a,b,alternative="greater")$p.value,error=function(e)1)
    data.frame(LS=x,lab=ord$lab[i],col=L,mean_pct=mean(a)*100,p=p)}))}))
enr$fdr<-p.adjust(enr$p,"BH"); enr$neglog10fdr<- -log10(enr$fdr)
enr$lab<-factor(enr$lab,levels=rev(lv)); enr$col<-factor(enr$col,levels=F3$CBF)
szmax<-ceiling(max(enr$mean_pct))
pc <- ggplot(enr,aes(col,lab,size=mean_pct,fill=neglog10fdr))+
  geom_hline(yintercept=yline,colour="grey40",linewidth=0.5)+geom_point(shape=21,colour="grey30",stroke=0.3)+
  scale_fill_seq(limits=c(0,2.5),oob=squish,name=expression(-log[10]~FDR))+
  scale_size(range=c(1,11),limits=c(0,szmax),name="mean % (deconv.)")+
  labs(x=NULL,y=NULL)+ggtitle("CBF cytogenetic")+theme_minimal(base_size=10)+
  theme(axis.text.x=element_text(angle=45,hjust=1,face="bold"),axis.text.y=element_text(size=9),
        panel.grid.major=element_line(colour="grey93"))
blk<-data.frame(y=c(nrow(ord)-nfav/2+0.5,(nrow(ord)-nfav)/2+0.5),lab=c("Favorable","Poor"))
pblk<-ggplot(blk,aes(1,y,label=lab,colour=lab))+geom_text(angle=90,fontface="bold",size=4)+
  scale_colour_manual(values=c(Favorable=unname(PROG["favorable"]),Poor=unname(PROG["poor"])),guide="none")+
  scale_y_continuous(limits=c(0.5,nrow(ord)+0.5),expand=expansion(mult=0.06))+
  scale_x_continuous(limits=c(0.5,1.5))+coord_cartesian(clip="off")+theme_void()
prs<-ord; prs$lab<-factor(prs$lab,levels=rev(lv)); prs$stem<-stemv[prs$LS]
ps<-ggplot(prs,aes("HSC/LSC",lab,fill=stem))+geom_hline(yintercept=yline,colour="grey40",linewidth=0.5)+
  geom_tile(colour="white",linewidth=0.4)+scale_fill_div(midpoint=0,limits=c(-1,2),oob=squish,name="stemness\n(z)")+
  labs(x=NULL,y=NULL)+theme_minimal(base_size=10)+
  theme(axis.text.y=element_blank(),axis.text.x=element_text(angle=45,hjust=1,size=9),panel.grid=element_blank())
prs$txt<-sprintf("%.2f  (p=%s)",prs$HR,formatC(prs$p,format="fg",digits=2)); prs$best<-prs$p==min(prs$p)
ptxt<-ggplot(prs,aes(1,lab,label=txt,fontface=ifelse(best,"bold","plain")))+
  geom_hline(yintercept=yline,colour="grey40",linewidth=0.5)+geom_text(size=3.1,colour="grey15")+
  scale_x_continuous(limits=c(0.4,1.6))+ggtitle("EFS Cox  HR (p)")+coord_cartesian(clip="off")+theme_void()+
  theme(plot.title=element_text(size=9,hjust=0.5))
fig<-pblk+pc+ps+ptxt+plot_layout(widths=c(0.05,0.5,0.12,0.34),guides="collect")+
  plot_annotation(theme=theme(plot.margin=margin(16,8,10,8)))
ggsave(file.path(F3$OUT,"panelA_cbf_states.pdf"),fig,width=9.5,height=6.6)
ggsave(file.path(F3$OUT,"panelA_cbf_states.png"),fig,width=9.5,height=6.6,dpi=200)
cat(sprintf("wrote panelA_cbf_states: %d states\n",nrow(ord)))
