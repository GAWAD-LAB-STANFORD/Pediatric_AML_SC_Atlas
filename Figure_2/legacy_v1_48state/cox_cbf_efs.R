suppressMessages(library(survival))
SL="/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad"
PKG="/Users/chuckgawad/Desktop/ALSF_AML_2026_updated"
d=read.csv(file.path(SL,"target_ls/LS48_survival.csv"))
ef=read.csv(file.path(PKG,"__SUBMISSION_PACKAGE/08_Source_Data/Figure_2/Figure2E_TARGET_survival_input.csv"))[,c("sample","efs_time","efs_event")]
d=merge(d,ef,by="sample")
cbf=d[d$subtype2 %in% c("t(8;21)","inv(16)") & !is.na(d$efs_time) & !is.na(d$efs_event),]
cat(sprintf("CBF with EFS: %d patients, EFS events=%d\n\n",nrow(cbf),sum(cbf$efs_event)))
ann=read.csv(file.path(SL,"ls48/LS48_annotation.csv"))[,c("LS","lineage_program")]
stem=read.csv(file.path(SL,"ls48/LS48_stem_scores.csv"))[,c("LS","stem_z")]
ls=paste0("LS_",1:48); res=data.frame()
for(x in ls){ fr=cbf[[x]]; if(sd(fr)==0||sum(fr>0)<10) next
  z=scale(fr)[,1]; f=tryCatch(coxph(Surv(efs_time,efs_event)~z,data=cbf),error=function(e)NULL); if(is.null(f))next
  s=summary(f); res=rbind(res,data.frame(LS=x,HR=round(s$conf.int[1],3),lo=round(s$conf.int[3],2),hi=round(s$conf.int[4],2),p=s$coefficients[5]))}
res$fdr=p.adjust(res$p,"BH"); res$dir=ifelse(res$HR>1,"poor","favorable")
res=merge(res,ann,by="LS"); res=merge(res,stem,by="LS"); res$LSc=gsub("_","",res$LS); res=res[order(res$p),]
options(width=160); cat("Per-state Cox on EFS within CBF (top):\n")
print(head(res[,c("LSc","dir","HR","lo","hi","p","fdr","lineage_program","stem_z")],10),row.names=FALSE)
cat(sprintf("\nEFS within CBF: FDR<0.10 = %d, nominal p<0.05 = %d\n",sum(res$fdr<0.10),sum(res$p<0.05)))
