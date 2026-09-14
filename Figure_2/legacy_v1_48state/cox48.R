suppressMessages(library(survival))
SL="/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad/target_ls"
d=read.csv(file.path(SL,"LS48_survival.csv")); d=d[!is.na(d$os_time)&!is.na(d$os_event),]
ls=paste0("LS_",1:48); res=data.frame()
for(x in ls){fr=d[[x]]; if(sd(fr)==0) next; z=scale(fr)[,1]; f=coxph(Surv(os_time,os_event)~z,data=d); s=summary(f)
 res=rbind(res,data.frame(LS=x,HR=round(s$conf.int[1],3),p=s$coefficients[5]))}
res$fdr=p.adjust(res$p,"BH"); res$dir=ifelse(res$HR>1,"poor","favorable"); res=res[order(res$p),]
write.csv(res,file.path(SL,"LS48_prognosis.csv"),row.names=FALSE)
cat(sprintf("prognosis-associated: FDR<0.05 = %d, FDR<0.10 = %d, nominal p<0.05 = %d (of 48)\n",sum(res$fdr<0.05),sum(res$fdr<0.10),sum(res$p<0.05)))
options(width=150); print(head(res[res$fdr<0.10,],20),row.names=FALSE)
