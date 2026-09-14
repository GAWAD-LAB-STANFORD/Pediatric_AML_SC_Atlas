# ====================================================================
# CIBERSORTx_Prognosis_risk_T_NK.R  |  CD96 figure pipeline component
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript CIBERSORTx_Prognosis_risk_T_NK.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

rm(list=ls())
suppressPackageStartupMessages({
  library(tidyverse)
})
library(ggplot2)
library(data.table)
library(dplyr)
library(stringr)
library(ggsignif)
library(MetBrewer)
library(GENIE3)
library("readxl")
library(biomaRt)
library(immunedeconv)
getwd()
generate_mat <- function(path, pattern) {
  files = list.files(path, pattern, full.names = TRUE, recursive=TRUE, include.dirs=TRUE)
  mat = as.data.frame(do.call(rbind, lapply(files, function(x) fread(x, stringsAsFactors = FALSE))))
  mat = as.data.frame(mat)
  return(mat)
}

NCI_CIBERSORT <- generate_mat("/Users/yakun/Downloads/CIBERSORTx/CIBERSORTx_T_NK", "\\_Results.csv$")

meta_1 <-read_excel("/Users/yakun/Downloads/ALSF_AML_plot/NCBI_PAML/TARGET_AML_ClinicalData_AAML03P1_AAML0531_CCG2961_20211201.xlsx")
meta_2 <-read_excel("/Users/yakun/Downloads/ALSF_AML_plot/NCBI_PAML/TARGET_AML_ClinicalData_AML1031_20211201.xlsx")
meta_3 <-read_excel("/Users/yakun/Downloads/ALSF_AML_plot/NCBI_PAML/TARGET_AML_ClinicalData_AAML1031_AAML0631_additionalCasesForSortedCellsAndCBExperiment_20220330.xlsx")

colnames(meta_1)
colnames(meta_2)
colnames(meta_3)
meta<-rbind(meta_1,meta_2,meta_3)
meta<-meta%>%separate(col="TARGET USI",sep = "-",c("Project","ProjectID","Patient_ID","Sample_Type"),extra="merge")

NCI_df<-NCI_CIBERSORT%>%separate(col="Mixture",sep = "-",c("Project","ProjectID","Patient_ID","Sample_Type"),extra="merge")

unique(NCI_df$Sample_Type)
#primary leukemia
NCI_df<-filter(NCI_df, Sample_Type%in%c("09A","03A",
                                        "Unsorted-09A"))
#Recurrent leukemia
#NCI_df<-filter(NCI_df, Sample_Type%in%c("04A",
 #                                       "40A"))

NCI_df_RemoveDup<-NCI_df%>% distinct(Patient_ID, .keep_all = TRUE)
NCI_df_RemoveDup<-unite(NCI_df_RemoveDup,PatientID,Project,ProjectID,Patient_ID, sep = "-", remove = TRUE, na.rm = FALSE)

#NCI_df_RemoveDup$PatientID<-sub("\\-", "\\.", NCI_df_RemoveDup$PatientID)
Sample_List<-NCI_df_RemoveDup$PatientID


NCI_df_LSC<- NCI_df_RemoveDup%>%mutate(LSC = rowSums(across(starts_with("LSC_")), na.rm=TRUE))

meta_f<-meta[,-4]
meta_f<-unite(meta_f,PatientID,Project,ProjectID,Patient_ID, sep = "-", remove = TRUE, na.rm = FALSE)
meta_f<-filter(meta_f,`Vital Status`%in%c("Dead","Alive"))
meta_f <-meta_f %>%mutate(Prognosis=0)%>%mutate(Prognosis=ifelse(`Vital Status`=='Alive',0,Prognosis))%>%
  mutate(Prognosis=ifelse(`Vital Status`=="Dead",1,Prognosis))
meta_f$`Bone marrow leukemic blast percentage (%)`<-as.numeric(as.character(meta_f$`Bone marrow leukemic blast percentage (%)`))
#meta_f<-filter(meta_f, `Bone marrow leukemic blast percentage (%)`>40)


colnames(meta_f)
Events_OS_df<-subset(meta_f,select=c("PatientID","Event Free Survival Time in Days",
                                     "Overall Survival Time in Days",
                                     "Bone marrow leukemic blast percentage (%)",
                                     "Year of Diagnosis","Year of Last Follow Up","Vital Status","First Event"))
Events_OS_df[,2:(ncol(Events_OS_df)-2)]<-lapply(Events_OS_df[,2:(ncol(Events_OS_df)-2)], function(x) as.numeric(as.character(x)))

Events_OS_df<-Events_OS_df%>%mutate(status=`Vital Status`)
Events_OS_df<-Events_OS_df%>%mutate(status=ifelse(`Vital Status`=='Alive', 'Censored', status))
Events_OS_df<-Events_OS_df%>%mutate(time=`Year of Last Follow Up`-`Year of Diagnosis`)


df_final<-merge(NCI_df_RemoveDup[,-2], Events_OS_df,by='PatientID', all=F)
colnames(df_final)

df_final<-dplyr::rename(df_final, `E.S.`=`Event Free Survival Time in Days`)
df_final<-dplyr::rename(df_final, `O.S.`=`Overall Survival Time in Days`)
df_final<-dplyr::rename(df_final, `Blast`=`Bone marrow leukemic blast percentage (%)`)


df_final<-unique(df_final)
df_final<-filter(df_final, !(time<3 & status=='Censored'))
data<-df_final[,-c(1,19,20,21,24:30)]


data<-data%>%mutate(NaïveCD4_sum=data$`Naïve CD4 T`+data$`AML-Naïve CD4 T`)
data<-data%>%mutate(NaïveCD8_sum=data$`Naïve CD8 T`+data$`AML-Naïve CD8 T`)
data<-data%>%mutate(Activated.CD4_sum=data$`T reg`+data$`AML-T reg`)
data<-data%>%mutate(GZMBCD8_sum=data$`GZMB DN T`+data$`GZMB CD8 T`+data$`AML-GZMB CD8 T`)
data<-data%>%mutate(GZMKCD8_sum=data$`GZMK CD8 T`+data$`GZMK CD8 T`+data$`AML-GZMK CD8 T`)
data<-data%>%mutate(CTL_sum=data$`GZMB DN T`+data$`GZMB CD8 T`+data$`AML-GZMB CD8 T`+data$`GZMK CD8 T`+data$`GZMK CD8 T`+data$`AML-GZMK CD8 T`)
data<-data%>%mutate(AML_CTL_sum=data$`AML-GZMB CD8 T`+data$`AML-GZMK CD8 T`)
data<-data%>%mutate(MAIT_sum=data$MAIT+data$`AML-MAIT`)
data<-data%>%mutate(GZMBNK_sum=data$`GZMB NK`+data$`AML-NK`)
data<-data%>%mutate(GZMKNK_sum=data$`GZMK NK`)
data[,1:(ncol(data))]<-lapply(data[,1:(ncol(data))], function(x) as.numeric(as.character(x)))
data[is.na(data)] <- 0.00

data<-data[,c(18:(ncol(data)))]
rownames(data)<-df_final$PatientID


library(corrplot)
testRes<-cor.mtest(data,conf.level = 0.95)
M = cor(data)

col=colorRampPalette(met.brewer("Signac", 2))(200)

corrplot(M, p.mat = testRes$p, method = 'circle',insig='blank',
         type ='upper',
         diag = FALSE,
         tl.cex = 0.8,addCoef.col = 1, tl.col = "brown",   # Change font size of correlation coefficients
         number.cex = 0.6,
         tl.srt = 45,
         order = 'alphabet',
)$corrPos -> p1

library(corrr)
Focus="O.S."
df_corr<-data %>% correlate() %>% focus(Focus)

good_corr<-filter(df_corr, O.S.>0)

poor_corr<-filter(df_corr, O.S.< 0)

paste(shQuote(good_corr$term), collapse=", ")

good_list<-c('GZMK NK', 'T reg', 
             'Naïve CD8 T', 'GZMB DN T', 
             'AML-GZMK CD8 T', 'GZMB CD8 T',
             'AML-T reg', 'AML-NK', 'AML-MAIT', 
             'E.S.', 'Activated.CD4_sum', 'GZMKCD8_sum', 'MAIT_sum', 'GZMKNK_sum')
poor_list<-c('GZMB NK', 'Naïve CD4 T',
             'GZMK CD8 T', 
             'MAIT', 'AML-Naïve CD4 T',
             'Effector memory CD8 T', 'AML-Naïve CD8 T',
             'AML-GZMB CD8 T', 'NaïveCD4_sum',
             'NaïveCD8_sum', 'GZMBCD8_sum', 
             'CTL_sum', 'GZMBNK_sum')



poor_list<-unique(poor_list)
keep_list<-poor_list

data_2<-subset(data,select = keep_list)

library("PerformanceAnalytics")
my_data <- data
chart.Correlation(my_data, histogram=TRUE, pch=19)

###############################################
library("survival")
library("survminer")
library(ranger)
library(ggfortify)
library(MetBrewer)

Events_OS_df<-subset(meta_f,select=c("PatientID","Event Free Survival Time in Days",
                                     "Overall Survival Time in Days",
                                     "Bone marrow leukemic blast percentage (%)",
                                     "Year of Diagnosis","Year of Last Follow Up","Vital Status","First Event"))
Events_OS_df[,2:(ncol(Events_OS_df)-2)]<-lapply(Events_OS_df[,2:(ncol(Events_OS_df)-2)], function(x) as.numeric(as.character(x)))

Events_OS_df<-Events_OS_df%>%mutate(status=`Vital Status`)
Events_OS_df<-Events_OS_df%>%mutate(status=ifelse(`Vital Status`=='Alive', 'Censored', status))
Events_OS_df<-Events_OS_df%>%mutate(time=`Year of Last Follow Up`-`Year of Diagnosis`)
Events_OS_df<-Events_OS_df%>%mutate(status= 2 )
Events_OS_df<-Events_OS_df%>%mutate(status=ifelse(`Vital Status`=='Alive', 1, status))


df_final<-merge(NCI_df_RemoveDup[,-2], Events_OS_df,by='PatientID', all=F)
colnames(df_final)

df_final<-dplyr::rename(df_final, `E.S.`=`Event Free Survival Time in Days`)
df_final<-dplyr::rename(df_final, `O.S.`=`Overall Survival Time in Days`)
df_final<-dplyr::rename(df_final, `Blast`=`Bone marrow leukemic blast percentage (%)`)
df_final<-filter(df_final, !(time<3 & `Vital Status`!='Dead'))
df_final<-unique(df_final)
data<-df_final

data<-data%>%mutate(NaïveCD4_sum=data$`Naïve CD4 T`+data$`AML-Naïve CD4 T`)
data<-data%>%mutate(NaïveCD8_sum=data$`Naïve CD8 T`+data$`AML-Naïve CD8 T`)
data<-data%>%mutate(Activated.CD4_sum=data$`T reg`+data$`AML-T reg`)
data<-data%>%mutate(GZMBCD8_sum=data$`GZMB DN T`+data$`GZMB CD8 T`+data$`AML-GZMB CD8 T`)
data<-data%>%mutate(GZMKCD8_sum=data$`GZMK CD8 T`+data$`AML-GZMK CD8 T`)
data<-data%>%mutate(CTL_sum=data$`GZMB DN T`+data$`GZMB CD8 T`+data$`GZMK CD8 T`)
data<-data%>%mutate(AML_CTL_sum=data$`AML-GZMB CD8 T`+data$`AML-GZMK CD8 T`)
data<-data%>%mutate(MAIT_sum=data$MAIT+data$`AML-MAIT`)
data<-data%>%mutate(GZMBNK_sum=data$`GZMB NK`+data$`AML-NK`)
data<-data%>%mutate(GZMKNK_sum=data$`GZMK NK`)

data[,-c(1,27:29)]<-lapply(data[,-c(1,27:29)], function(x) as.numeric(as.character(x)))
data[is.na(data)] <- 0.00

hist(data$CTL_sum)
quantile(data$CTL_sum,
         na.rm = T,
         probs = c(0.95)
         )
mean(data$CTL_sum)
colnames(data)

data$Group<-'CTL Low'
data<-data%>%mutate(Group=ifelse(CTL_sum >0.08,
                                    'CTL High', Group))

fit <- survfit(Surv(`O.S.`,status) ~ Group, data =data)

print(fit)
# Summary of survival curves
summary(fit)
# Access to the sort summary table
summary(fit)$table

d <- data.frame(time = fit$time,
                n.risk = fit$n.risk,
                n.event = fit$n.event,
                n.censor = fit$n.censor,
                surv = fit$surv,
                upper = fit$upper,
                lower = fit$lower
)


head(d)
ggsurvplot(fit,
           data =data,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
#####################################################
cluster_list<-c("`GZMB NK`" ,              
                "`GZMK NK`",                "`T reg`" ,                
                "`Naïve CD4 T`" ,           "`Naïve CD8 T`",           
                "`GZMK CD8 T`"  ,           "`GZMB DN T`" ,            
                "MAIT"    ,               "`AML-Naïve CD4 T`" ,      
                "`AML-GZMK CD8 T`",         "`GZMB CD8 T`" ,           
                "`Effector memory CD8 T`",  "`AML-Naïve CD8 T`",       
                "`AML-T reg`",              "`AML-NK`" ,               
                "`AML-GZMB CD8 T`" ,        "`AML-MAIT`",
                 "NaïveCD4_sum",
                "NaïveCD8_sum", "Activated.CD4_sum",
                "GZMBCD8_sum",   "GZMKCD8_sum", "AML_CTL_sum", 
                "CTL_sum",      "MAIT_sum",
                "GZMBNK_sum",      "GZMKNK_sum")

cluster_list<-c(           
              "`Effector memory CD8 T`", 
                "NaïveCD4_sum",
                "NaïveCD8_sum", "Activated.CD4_sum",
                "GZMBCD8_sum",   "GZMKCD8_sum", "AML_CTL_sum", 
                "CTL_sum",      "MAIT_sum",
                "GZMBNK_sum",      "GZMKNK_sum")
res.cox <- coxph(Surv(time, status) ~ CTL_sum, data = data)
res.cox
summary(res.cox)
#假设我们要对如下5个特征做单因素cox回归分析
covariates <- cluster_list
#分别对每一个变量，构建生存分析的公式

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
#循环对每一个特征做cox回归分析
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data)})
#提取HR，95%置信区间和p值
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         #获取HR
                         HR <-signif(x$coef[2], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })
#转换成数据框，并转置
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)

colnames(res)

res<-res%>%separate(col="HR (95% CI for HR)", sep=' ',c('HR','95% CI for HR'))
res[,1:2]<-lapply(res[,1:2], function(x) as.numeric(as.character(x)))
res<-filter(res, p.value<0.005)
res_HR<-filter(res, HR>1)
rownames(res_HR)
res_LR<-filter(res, HR<1)
rownames(res_LR)
write.table(file="T_NK_univariate_cox_result_primary_leukemia.txt",as.data.frame(res),quote=F,sep="\t")
###################################################
library(ggrisk)
library(rms)
colnames(data)
data[,31:39]<-lapply(data[,31:39], function(x) x*100)

data<-data%>%mutate(status= 1 )

data<-data%>%mutate(status=ifelse(`Vital Status`=='Alive', 0, status))
colnames(data)
fit <- cph(Surv(O.S.,status)~MAIT_sum+CTL_sum,
           data[,c(23,29,36,37)])

ggrisk(fit,
       code.highrisk = 'High Risk',#高风险标签，默认为 ’High’
       code.lowrisk = 'Low Risk', #低风险标签，默认为 ’Low’
       title.A.ylab='Risk Score', #A图 y轴名称
       title.B.ylab='Survival Time(Days)', #B图 y轴名称，注意区分year month day
       title.A.legend='Risk Group', #A图图例名称
       title.B.legend='Status',     #B图图例名称
       title.C.legend='Fractions', #C图图例名称
       relative_heights=c(0.1,0.1,0.01,0.15), #A、B、热图注释和热图C的相对高度    
       color.A=c(low='#023FA5',high='#FF7F0E'),#A图中点的颜色
       color.B=c(code.0='#FF7F0E',code.1='#023FA5'), #B图中点的颜色
       color.C=c(low='#023FA5',median='white',high='#FF7F0E'), #C图中热图颜色
       vjust.A.ylab=1, #A图中y轴标签到y坐标轴的距离,默认是1
       vjust.B.ylab=2  #B图中y轴标签到y坐标轴的距离,默认是2)
)


death_times <- r_fit$unique.death.times 
surv_prob <- data.frame(r_fit$survival)
avg_prob <- sapply(surv_prob,mean)