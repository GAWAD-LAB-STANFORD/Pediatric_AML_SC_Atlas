# ====================================================================
# FigureS3.LSC_fraction_in_leiden.R  |  Figure S4 (LSC fraction; risk model)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : ALSF_AML_leiden.csv, AML-sample-infor.csv, Gene_feature_score_Combo_enriched.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript FigureS3.LSC_fraction_in_leiden.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
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
library(ggrepel)


setwd("~/Downloads/ALSF_AML_plot")
df_GF <- read.delim2("Gene_feature_score_Combo_enriched.csv", sep = ",", header = TRUE)
colnames(df_GF)
df_GF<- df_GF%>%mutate(across(c(HSC_Score:LSC_EPPERT_Score), function(x) as.numeric(as.character(x))))

df_LSC<- filter(df_GF,LSC_EPPERT_Score>0.2)

df_GF_sum <- df_GF %>% 
  group_by(leiden) %>% 
  summarise(leiden_total = n())

df_LSC_sum <- df_LSC %>% 
  group_by(leiden) %>% 
  summarise(leiden_LSC = n())

df<-merge(df_GF_sum,df_LSC_sum, by='leiden',all = TRUE)

df_P<-df%>%mutate(Percentage = leiden_LSC/leiden_total)

rownames(df_P)<-df_P$leiden

ggplot(df_P, aes(x=leiden, y=Percentage)) +
  geom_point()+
  geom_text_repel(label=rownames(df_P))+
  theme_minimal()+
  ylab("LSC_EPPERT_Score>0.2 Fraction")+xlab("Leiden")+
  geom_hline(yintercept=0.1, linetype="dashed", color = "red")+
  theme(plot.title = ggplot2::element_text(hjust = 0.5))+
  theme(axis.text.x =element_text(face="bold", 
                                  size=10),
        axis.text.y = element_text(face="bold", 
                                   size=10))

df_P_top<-filter(df_P, Percentage>0.1)


################################################LSC17>0.0#########################
df_LSC<- filter(df_GF,LSC17_Score>0.0)

df_LSC_sum <- df_LSC %>% 
  group_by(leiden) %>% 
  summarise(leiden_LSC = n())

df<-merge(df_GF_sum,df_LSC_sum, by='leiden',all = TRUE)

df_P_LSC17<-df%>%mutate(Percentage = leiden_LSC/leiden_total)

rownames(df_P_LSC17)<-df_P_LSC17$leiden

ggplot(df_P_LSC17, aes(x=leiden, y=Percentage)) +
  geom_point()+
  geom_text_repel(label=rownames(df_P_LSC17))+
  theme_minimal()+
  ylab("LSC17_Score>0.0 Fraction")+xlab("Leiden")+
  geom_hline(yintercept=0.1, linetype="dashed", color = "red")+
  theme(plot.title = ggplot2::element_text(hjust = 0.5))+
  theme(axis.text.x =element_text(face="bold", 
                                  size=10),
        axis.text.y = element_text(face="bold", 
                                   size=10))

df_P_LSC17<-filter(df_P_LSC17, Percentage > 0.1)

############################################
setwd("~/Downloads/ALSF_AML_plot")
df_clic <- read.delim2("AML-sample-infor.csv", sep = ",", header = TRUE)
df_clic <-df_clic %>%mutate(Prognosis=0)%>%mutate(Prognosis=ifelse(Alive.deceased=='Alive',0,Prognosis))%>%
  mutate(Prognosis=ifelse(Alive.deceased=='Deceased',1,Prognosis))

df_clic$MRD<-df_clic$MRD...
df_clic$Blast<-df_clic$Blast..
df_clic<-subset(df_clic,select=c(2,36:44))
colnames(df_clic)
df_clic <- df_clic%>%mutate(across(c(Alive:Blast), function(x) as.numeric(as.character(x))))
df_clic <- df_clic%>%mutate(across(c(Alive:Blast), ~replace(., is.na(.), 0)))

df_clic <- filter(df_clic, !Sample%in%c('AML4363','AML647'))

df_GF <- read.delim2("ALSF_AML_leiden.csv", sep = ",", header = TRUE)
colnames(df_GF)
df_GF<- df_GF%>%mutate(across(c(X1:X101), function(x) as.numeric(as.character(x))))

df_GF <- df_GF%>%mutate(across(starts_with("X"), ~ ifelse( .x>0, .x/df_GF$sum, 0), .names="Fraction_{col}"))

df_GF <- df_GF%>%mutate(across(starts_with("X"), ~ ifelse( .x>0, .x/df_GF$sum, 0), .names="Fraction_{col}"))
df_GF<-df_GF %>% separate(col="Sample", sep = "-", c("Sample", "Drop1"))

df<-merge(df_clic, df_GF[,-2], by='Sample')

df[,2:(ncol(df))]<-lapply(df[,2:(ncol(df))], function(x) as.numeric(as.character(x)))

data<-df[,-1]
rownames(data)<-df[,1]
data =data[, !(colSums(data) == 0)]
data =data[, !(colSums(data) == 26)]
library(mccf1)

test_list <- c(colnames(data[, 113:(ncol(data))]))
test_list 

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$Prognosis, data[ , t])
  df_mccf <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df_mccf)
}

colnames(mccf_final)

factor_list<-dplyr::filter(mccf_final, normalized_mcc>0.6)
factor_list<-dplyr::filter(factor_list, f1>0.55)
factor_list<-unique(factor_list$Factor)
paste(shQuote(factor_list), collapse=", ")

