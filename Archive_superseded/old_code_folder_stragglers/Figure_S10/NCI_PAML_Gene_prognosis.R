# ====================================================================
# NCI_PAML_Gene_prognosis.R  |  Figure S10 (surfaceome MCC/F1 screen)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : NCI_TARGET_gene_counts_filtered_anno.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript NCI_PAML_Gene_prognosis.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

rm(list=ls())
suppressPackageStartupMessages({
  library(tidyverse)
})
library("readxl")
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(dplyr)
library(tibble)
library(pROC)
library(mccf1)
library(MetBrewer)
setwd("~/Downloads/ALSF_AML_plot/NCBI_PAML")
####################################################################################################
df_NCBI_F<-read.delim2("NCI_TARGET_gene_counts_filtered_anno.csv", header=T, sep=',')
NCI_df<-df_NCBI_F

NCI_df$PatientID_SampleType<-df_NCBI_F$Sample
NCI_df<-NCI_df%>%separate(col="PatientID_SampleType",sep = "\\.",c("Project","ProjectID","Patient_ID","Sample_Type"),extra="merge")
NCI_df$Sample_Type
NCI_df_Primary<-filter(NCI_df, Sample_Type%in%c("09A",
                                                "03A",
                                                "Sorted-leukemic-09A",
                                                "Unsorted-09A"))

NCI_df_Primary<-unite(NCI_df_Primary,
                      PatientID,
                      Project,
                      ProjectID,
                      Patient_ID, 
                      sep = "-", 
                      remove = FALSE, na.rm = FALSE)


meta_1 <-read_excel("/Users/yakun/Downloads/ALSF_AML_plot/NCBI_PAML/TARGET_AML_ClinicalData_AAML03P1_AAML0531_CCG2961_20211201.xlsx")
meta_2 <-read_excel("/Users/yakun/Downloads/ALSF_AML_plot/NCBI_PAML/TARGET_AML_ClinicalData_AML1031_20211201.xlsx")
meta_3 <-read_excel("/Users/yakun/Downloads/ALSF_AML_plot/NCBI_PAML/TARGET_AML_ClinicalData_AAML1031_AAML0631_additionalCasesForSortedCellsAndCBExperiment_20220330.xlsx")

colnames(meta_1)
colnames(meta_2)
colnames(meta_3)
meta<-rbind(meta_1,meta_2,meta_3)
meta<-meta%>% distinct(`TARGET USI`, .keep_all = TRUE)
names(meta)[names(meta) == 'TARGET USI'] <- 'PatientID'


meta_f<-meta[,-4]
meta_f$`Bone marrow leukemic blast percentage (%)`<-as.numeric(as.character(meta_f$`Bone marrow leukemic blast percentage (%)`))
#meta_f<-filter(meta_f, `Bone marrow leukemic blast percentage (%)`>40)
meta_f_RemoveDup<-meta_f%>% distinct(PatientID, .keep_all = TRUE)

colnames(meta_f)
###############################
library(edgeR)
library(ggplot2)
NCI_df_Primary<-NCI_df_Primary%>% distinct(PatientID, .keep_all = TRUE)
df_NCI_F<-t(NCI_df_Primary[, -c((ncol(NCI_df_Primary)-6):ncol(NCI_df_Primary))])

rownames(df_NCI_F)<-colnames(NCI_df_Primary[, -c((ncol(NCI_df_Primary)-6):ncol(NCI_df_Primary))])

colnames(df_NCI_F)<-NCI_df_Primary$PatientID
df_NCI_F<-as.data.frame(df_NCI_F)
df_NCI_F[,1:ncol(df_NCI_F)]<-lapply(df_NCI_F[,1:ncol(df_NCI_F)],function(x) as.numeric(as.character(x)))

count_matrix <-as.matrix(df_NCI_F)

dge <- DGEList(counts = count_matrix, group = factor(colnames(df_NCI_F)))

keep <- filterByExpr(y = dge)

dge <- dge[keep, , keep.lib.sizes=FALSE]

dge <- calcNormFactors(object = dge)

cpms <- as.data.frame(cpm(dge, log=FALSE)) 

######################################Cox Regression#######
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

Events_OS_df<-dplyr::rename(Events_OS_df, `E.S.`=`Event Free Survival Time in Days`)
Events_OS_df<-dplyr::rename(Events_OS_df, `O.S.`=`Overall Survival Time in Days`)
Events_OS_df<-dplyr::rename(Events_OS_df, `Blast`=`Bone marrow leukemic blast percentage (%)`)

Events_OS_df<-filter(Events_OS_df, !(`First Event`=='Censored' & `O.S.`<200))

Events_OS_df<-as.data.frame(Events_OS_df)
rownames(Events_OS_df) <- 1:nrow(Events_OS_df)
colnames(Events_OS_df)


#####################
library("survival")
library("survminer")
library(ranger)
library(ggfortify)
library(MetBrewer)
library(cowplot)

df_pro<-as.data.frame(t(cpms))

rownames(df_pro)<-colnames(cpms)
colnames(df_pro)<-rownames(cpms)
df_pro$PatientID<-rownames(df_pro)


df_pro_input<-merge(Events_OS_df,df_pro,by='PatientID')

CSF<-c('PTAFR', 'COL9A2', 'ARTN', 'FCGR1A', 'ECM1', 'BGLAP',
       'FCGR2A', 'IL5RA', 'LAMB2', 'STAB1', 'CD96', 'SLC15A2',
       'TLR2', 'CD180', 'ADRB2', 'CSF1R', 'HRH2', 'F13A1', 
       'LY86', 'THSD7A', 'ANGPT2', 'MAMDC2', 'PTGDS', 'NRXN2',
       'RNASE3', 'SERPINA1', 'AMN', 'LTK', 'NRG4', 'MSLN', 
       'TPSAB1', 'TPSD1', 'ZG16B', 'ITGAX', 'CES1', 'ICAM1', 
       'APOC2', 'SIGLEC12', 'VSTM1', 'OSCAR', 'LILRA5', 
       'LILRB4', 'JAG1', 'CD93', 'LAMA5', 'GGT5', 'CSF2RA',
       'AR', 'IL13RA1')

test_list<-CSF   # (moved below the CSF definition; was used before assignment)

test_list<-CSF

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(df_pro_input$status, df_pro_input[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}


write_csv(mccf_final, 'MCCF1_AML_specific_gene_list_NCI_PAML_validation.csv')

mccf_plot<-filter(mccf_final, normalized_mcc>0.54)

List<-unique(mccf_plot$Factor)
List

paste(shQuote(List), collapse=", ")

intersect(List,AML_module_list$gene_name)


df_pro_input$Group<-'Low'
df_pro_input<-df_pro_input%>%
  mutate(Group=ifelse(SIGLEC12>17, 'High', Group))

fit <- survfit(Surv(`O.S.`,status) ~Group, data = df_pro_input)

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

p<-ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata",# Change risk table color by groups
           risk.table.fontsize=3,
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_gray(),# Change ggplot2 theme
           palette = "Dark2")
print(p)
p$plot <- p$plot + labs(title = "SIGLEC12 threashold=17(CPM)")
print(p)
#####################################################
res.cox <- coxph(Surv(O.S., status) ~ GGT5, data = df_pro_input )
res.cox
summary(res.cox)
#假设我们要对如下5个特征做单因素cox回归分析
covariates <- colnames(df_pro)
#分别对每一个变量，构建生存分析的公式

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(O.S., status)~', x)))
#循环对每一个特征做cox回归分析
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = df_pro_input)})
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
res_f<-filter(res, p.value<0.05)
res_HR<-filter(res_f, HR>1)
rownames(res_HR)
res_LR<-filter(res_f, HR<1)
rownames(res_LR)
write.table(file="univariate_cox_result.txt",as.data.frame(res),quote=F,sep="\t")
###################################################
library(ggrisk)
library(rms)
colnames(data)
data[,2:30]<-lapply(data[,2:30], function(x) x*100)

data<-data%>%mutate(status= 1 )

data<-data%>%mutate(status=ifelse(`Vital Status`=='Alive', 0, status))

fit <- cph(Surv(O.S.,status)~AML_27+AML_24+AML_7+AML_4+AML_3,
           data[,c(35,41,5,2,15,3,26)])

fit <- cph(Surv(O.S.,status)~AML_2+AML_11+AML_14+AML_17+AML_18+AML_23+AML_26,
           data[,c(35,41,10,14,9,23,20,22,27)])

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





