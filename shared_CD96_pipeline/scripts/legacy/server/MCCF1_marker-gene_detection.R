# ====================================================================
# MCCF1_marker-gene_detection.R  |  Figure S10 (MCC/F1 marker detection)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : ALSF_RNA-seq_Combo_obs.csv, ALSF_RNA-seq_expr_AML_0.2.csv, ALSF_AML_ALL_cell_LSC_Module_MEs_pct_matrix.csv, ALSF_AML_Combo_obs.csv, ALSF_RNA-seq_Combo_obs.csv, ALSF_RNA-seq_expr_AML_0.2.csv ...
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript MCCF1_marker-gene_detection.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

rm(list=ls())
suppressPackageStartupMessages({
  library(tidyverse)
})

library(ggplot2)
library(dplyr)
library(tibble)
library(ggsci)
library(scales)
library(MetBrewer)
library(caTools)
library(MLmetrics)

library(pROC)
library(mccf1)


HBM<-c('HealthyBM1','HealthyBM2','0_HealthyBM1','0_HealthyBM2')
######################################################################################
df_GE <- read.delim2("", sep = ",", header = TRUE)

df_GE[,2:(ncol(df_GE))]<-lapply(df_GE[,2:(ncol(df_GE))], as.numeric)
df_obs <- read.delim2("ALSF_RNA-seq_Combo_obs.csv", sep = ",", header = TRUE)
colnames(df_GE)
df_GE <- df_GE %>% mutate(across(c(TNFRSF4:IL13RA1), ~ ifelse(.>0,1,0)))


colnames(df_obs)

df_lineage<-subset(df_obs, select = c("X","lineage"))
df_lineage<-df_lineage%>%mutate(Value=1)
df_lineage<-spread(df_lineage, key="lineage", value = "Value")
colnames(df_lineage)
df_lineage <- df_lineage%>%mutate(across(c(`0_HSPC`:`T`), ~replace(., is.na(.), 0)))

df<-merge(df_GE, df_lineage,by="X")


data<-df[,-1]
rownames(data)<-df[,1]

test_list <- c(colnames(data[, 1:(ncol(data))]))
test_list 

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$AML, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}

factor_list<-c("COL23A1","SIGLEC12","FCGR1A",
               "NT5E",  "NTRK1",
               "TRGC2","TRGV9",
               "CAV1",
               "ADGRA2",
               "IL15RA",
               "MSLN",
               "TNFRSF4")
library(corrr)
df_corr<-data %>% correlate() %>% focus(AML)
colnames(mccf_final)
factor_list<-filter(mccf_final, normalized_mcc>0.55)
factor_list<-filter(factor_list, f1>0.25)
factor_list<-unique(factor_list$Factor)

paste(shQuote(factor_list), collapse=", ")

paste(factor_list, collapse=", ")

mccf_plot<-filter(mccf_final, Factor%in%factor_list)
p<-ggplot(mccf_plot, aes(x = f1, y = normalized_mcc, fill=Factor))+
  geom_point(size = 2, shape = 21)+ylim(0,1)+xlim(0,1)+
  theme_minimal()+
  geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
  theme(plot.title = ggplot2::element_text(hjust = 0.5))+
  ylab("Normalized MCC")+xlab("F1 Score")+
  theme(axis.text.x =element_text(face="bold", 
                                  size=12),
        axis.text.y = element_text(face="bold", 
                                   size=12))+
  coord_equal(ratio = 1)

p
###########################################################################
df_GE <- read.delim2("ALSF_RNA-seq_expr_HSPC_absent_0.01_AML_0.1.csv", sep = ",", header = TRUE)

df_GE[,2:(ncol(df_GE))]<-lapply(df_GE[,2:(ncol(df_GE))], as.numeric)
df_obs <- read.delim2("ALSF_RNA-seq_Combo_obs.csv", sep = ",", header = TRUE)
colnames(df_GE)
df_GE <- df_GE %>% mutate(across(c(TNFRSF4:IL13RA1), ~ ifelse(.>0,1,0)))


colnames(df_obs)

df_lineage<-subset(df_obs, select = c("X","Prognosis"))
df_lineage<-filter(df_lineage, !Prognosis=="0_HealthyBM")
df_lineage<-df_lineage%>%mutate(Prognosis_bi=1)%>%
  mutate(Prognosis_bi=ifelse(Prognosis=="Alive",0,Prognosis_bi))
colnames(df_lineage)
df<-merge(df_GE, df_lineage[,c(1,3)],by="X")

df_lineage<-subset(df_obs, select = c("X","lineage"))
df_lineage<-df_lineage%>%mutate(Value=1)
df_lineage<-spread(df_lineage, key="lineage", value = "Value")
df_lineage <- df_lineage%>%mutate(across(c(`0_HSPC`:`T`), ~replace(., is.na(.), 0)))

df<-merge(df_GE, df_lineage,by="X")
data<-df[,-1]
rownames(data)<-df[,1]

test_list <- c(colnames(data[, 1:(ncol(data))]))
test_list 

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$AML, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}

factor_list<-c("COL23A1","SIGLEC12","FCGR1A",
               "NT5E",  "NTRK1",
               "TRGC2","TRGV9",
               "CAV1",
               "ADGRA2",
               "IL15RA",
               "MSLN",
               "TNFRSF4")
library(corrr)
df_corr<-data %>% correlate() %>% focus(AML)
colnames(mccf_final)
factor_list<-filter(mccf_final, normalized_mcc>0.55)
factor_list<-filter(factor_list, f1>0.25)
factor_list<-unique(factor_list$Factor)

paste(shQuote(factor_list), collapse=", ")

paste(factor_list, collapse=", ")

mccf_plot<-filter(mccf_final, Factor%in%factor_list)
p<-ggplot(mccf_plot, aes(x = f1, y = normalized_mcc, fill=Factor))+
  geom_point(size = 2, shape = 21)+ylim(0,1)+xlim(0,1)+
  theme_minimal()+
  geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
  theme(plot.title = ggplot2::element_text(hjust = 0.5))+
  ylab("Normalized MCC")+xlab("F1 Score")+
  theme(axis.text.x =element_text(face="bold", 
                                  size=12),
        axis.text.y = element_text(face="bold", 
                                   size=12))+
  coord_equal(ratio = 1)

p
#'FCGR1A', 'KYNU', 'CDKN2A', 'MS4A6A', 'CLECL1', 'CCNA1', 'MSLN', 'PCAT18', 'APOC2', 'Prognosis_bi'
###############################################################################
df_GE <- read.delim2("ALSF_RNA-seq_expr_HSPC_Absent_0.01_0.1.csv", sep = ",", header = TRUE)
df_GE[,2:(ncol(df_GE))]<-lapply(df_GE[,2:(ncol(df_GE))], as.numeric)
df_obs <- read.delim2("ALSF_RNA-seq_Combo_obs.csv", sep = ",", header = TRUE)
colnames(df_GE)
df_GE <- df_GE %>% mutate(across(c(HES4:TBL1Y), ~ ifelse(.>0,1,0)))
#write.csv(df_GE, "ALSF_RNA-seq_expr_HSPC_Absent_binary.csv")

colnames(df_obs)

df_lineage<-subset(df_obs, select = c("X","Prognosis"))
df_lineage<-filter(df_lineage, !Prognosis=="0_HealthyBM")
df_lineage<-df_lineage%>%mutate(Prognosis_bi=1)%>%
  mutate(Prognosis_bi=ifelse(Prognosis=="Alive",0,Prognosis_bi))
colnames(df_lineage)
df<-merge(df_GE, df_lineage[,c(1,3)],by="X")


data<-df[,-1]
rownames(data)<-df[,1]

test_list <- c(colnames(data[, 1:(ncol(data))]))
test_list 

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$Prognosis_bi, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}


write.csv(mccf_final, "ALSF_RNA-seq_expr_AML_0.2_vs_Prognosis.csv")

library(corrr)
df_corr<-data %>% correlate() %>% focus(Prognosis_bi)
colnames(mccf_final)
factor_list<-filter(mccf_final, normalized_mcc>0.5)
factor_list<-filter(factor_list, f1>0.39)
factor_list<-unique(factor_list$Factor)
mccf_plot<-filter(mccf_final, Factor%in%factor_list)


p<-ggplot(mccf_plot, aes(x = f1, y = normalized_mcc, fill=Factor))+
  geom_point(size = 2, shape = 21)+ylim(0,1)+xlim(0,1)+
  theme_minimal()+
  geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
  theme(plot.title = ggplot2::element_text(hjust = 0.5))+
  ylab("Normalized MCC")+xlab("F1 Score")+
  coord_equal(ratio = 1)

p

paste(shQuote(factor_list), collapse=", ")
c('DHRS3', 'ARTN', 'FCGR1A', 'ECM1', 'S100A16', 'BGLAP', 'FCGR2A', 'SIPA1L2', 'AL391832.2', 
  'KYNU', 'PTH2R', 'SNORC', 'CD86', 'NDST3', 'IRX1', 'CD180', 'PITX1', 'COL23A1', 'FOXC1',
  'MDFI', 'NT5E', 'TRGC2', 'TRGV9', 'CAV1', 'TRBV28', 'TMEM176B', 'CDKN2A', 'AL392086.3', 
  'RASSF4', 'WT1', 'WT1.AS', 'PTPRJ', 'MS4A6A', 'MS4A4E', 'AP005273.2', 'PIWIL4', 'CLECL1', 
  'AC020656.1', 'CSRP2', 'CCNA1', 'POU4F1', 'MSLN', 'ITGAX', 'NETO2', 'IRX3', 'PPP1R27', 
  'PCAT18', 'TNFSF9', 'APOC2', 'SIGLEC12', 'LINC00649', 'MYO18B', 'GAS2L1', 'LGALS2', 'TSIX',
  'IL13RA1', 'Prognosis_bi') #mcc 0.55

c('FCGR1A', 'KYNU', 'CDKN2A', 'MS4A6A',
  'CLECL1', 'CCNA1', 'MSLN', 'PCAT18', 'APOC2', 'Prognosis_bi')#mcc 0.5 &f10.2
paste(factor_list, collapse=", ")
#########################################################
df_GE <- read.delim2("ALSF_RNA-seq_expr_HSPC.csv", sep = ",", header = TRUE)
df_obs <- read.delim2("ALSF_RNA-seq_Combo_obs.csv", sep = ",", header = TRUE)
colnames(df_GE)
df_GE <- df_GE %>% mutate(across(c(MFAP2:LDOC1), ~ ifelse(.>0,1,0)))
write.csv(df_GE, "ALSF_RNA-seq_expr_HSPC_Absent_binary.csv")

colnames(df_lineage)

df_lineage<-subset(df_obs, select = c("X","LSC_Cell_Type_2"))
df_lineage<-df_lineage%>%mutate(Value=1)
df_lineage<-spread(df_lineage, key="LSC_Cell_Type_2", value = "Value")
df_lineage <- df_lineage%>%mutate(across(c(`0_HSPC`:PlasmaB), ~replace(., is.na(.), 0)))

df<-merge(df_GE, df_lineage,by="X")


data<-df[,-1]
rownames(data)<-df[,1]

test_list <- c(colnames(data[, 1:(ncol(data))]))
test_list 

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$`0_HSPC`, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}

colnames(mccf_final)
factor_list<-filter(mccf_final, normalized_mcc>0.55)
factor_list<-filter(factor_list, f1>0.4)
factor_list<-unique(factor_list$Factor)

paste(shQuote(factor_list), collapse=", ")

paste(factor_list, collapse=", ")
#############################################################
df_GE <- read.delim2("ALSF_RNA-seq_expr_AML_0.2.csv", sep = ",", header = TRUE)
df_obs <- read.delim2("ALSF_RNA-seq_Combo_obs.csv", sep = ",", header = TRUE)
colnames(df_GE)
df_GE <- df_GE %>% mutate(across(c(NOC2L:VBP1), ~ ifelse(.>0,1,0)))
#write.csv(df_GE, "ALSF_RNA-seq_expr_HSPC_Absent_binary.csv")

colnames(df_obs)

df_lineage<-subset(df_obs, select = c("X","lineage"))
df_lineage<-df_lineage%>%mutate(Value=1)
df_lineage<-spread(df_lineage, key="lineage", value = "Value")
df_lineage <- df_lineage%>%mutate(across(c(`0_HSPC`:`T`), ~replace(., is.na(.), 0)))

df<-merge(df_GE, df_lineage,by="X")


data<-df[,-1]
rownames(data)<-df[,1]

test_list <- c(colnames(data[, 1:(ncol(data))]))
test_list 

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$`T`, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}
# "SH2D1B, FGFBP2, KLRF1, KLRD1, IL2RB, NK"
#"FCRLA, MME, DTX1, IGHV5-78, IGLC1, B"
# "SDC1, BHLHA15, IGHA2, TNFRSF17, TNFRSF13B, PlasmaB"
#"TGFBI, CD14, PLBD1, C5AR1, FPR1, Monocyte"
#"'PKLR', 'FAM178B', 'NMU', 'NECAB1', 'CDH1', 'Erythrocytes'"
#"CD8A, CD3G, CD27, PCED1B.AS1, GZMM, T"


write.csv(mccf_final, "ALSF_AML_MCCF_AML_AML_0.2.csv")
colnames(mccf_final)
factor_list<-filter(mccf_final, normalized_mcc>0.55)
factor_list<-filter(factor_list, f1>0.35)
factor_list<-unique(factor_list$Factor)
paste(factor_list, collapse=", ")
#"SH2D1B, KLRF1, KLRD1, KLRC1, S1PR5, AML-NK"
factor_list<-c("ENO1",
               "ITGA4",
               "CD74",
               "EGFL7",
               "CD44",
               "CD63",
               "SPN",
               "TYROBP",
               "MSN",
               'CSF3R')
factor_list<-c("TNFRSF14","CD34","ITGA6","TFPI", "CD38","PROM1","SLC38A1","ANPEP", "MBP","NECTIN2")
mccf_plot<-filter(mccf_final, Factor%in%factor_list)


p<-ggplot(mccf_plot, aes(x = f1, y = normalized_mcc, fill=Factor))+
  geom_point(size = 2.5, shape = 21)+ylim(0,1)+xlim(0,1)+
  theme_minimal()+
  geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
  theme(plot.title = ggplot2::element_text(hjust = 0.5))+
  ylab("Normalized MCC")+xlab("F1 Score")+
  theme(axis.text.x =element_text(face="bold", 
                                  size=10),
        axis.text.y = element_text(face="bold", 
                                   size=10))+
  coord_equal(ratio = 1)

p

paste(shQuote(factor_list), collapse=", ")



cellsurface<-read.table("GOCC_CELL_SURFACE.v2023.1.Hs.gmt", header=FALSE)
cellsurface_list<-as.list(cellsurface[1,c(3:931)])
intersect(factor_list,cellsurface_list)


###############################################################################
df_obs <- read.delim2("ALSF_RNA-seq_Combo_obs.csv", sep = ",", header = TRUE)
colnames(df_obs)
df_antibody<-subset(df_obs, select = c("X","lineage",
                                       "CD19_pos","CD274_pos","CD3_pos" ,"CD33_pos",                    
                                       "CD90_pos","CD10_pos","CD45RA_pos",                  
                                       "CD123_pos","CD7_pos" ,"CD49f_pos",                   
                                       "CD25_pos","CD279_pos","CD32_pos",                  
                                       "CD152_pos", "CD366_pos","CD235ab_pos",                
                                       "CD127_pos","CD71_pos","CD36_pos",                   
                                       "CD133_pos"))
df_antibody[,3:(ncol(df_antibody))]<-lapply(df_antibody[,3:(ncol(df_antibody))], as.numeric)

df_antibody_total_count<-df_antibody %>%
  group_by(SampleID) %>%
  summarize(Total = n())

df_antibody_pos <- df_antibody %>% group_by(SampleID) %>% 
  summarise(across(c(CD19_pos:Relapse_Score),sum),
            .groups = 'drop') %>%
  as.data.frame()

df_antibody_final <- merge(df_antibody_pos,df_antibody_total_count, by='SampleID')
df_antibody_final <-df_antibody_final %>%mutate(across(c(CD19_pos:Relapse_Score), function(x) x/df_antibody_final$Total))
df_antibody_final$Sample<-df_antibody_final$SampleID
data<-merge(df_antibody_final, df_clic, by='Sample')



rownames(data)<-data$Sample
data<-df_antibody[,-(1:2)]
rownames(data)<-df_antibody$X

test_list <- c(colnames(data[, 1:(ncol(data))]))
test_list 

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$CD49f_pos, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}

colnames(mccf_final)
factor_list<-filter(mccf_final, normalized_mcc>0.5)
factor_list<-unique(factor_list$Factor)

paste(shQuote(factor_list), collapse=", ")

paste(factor_list, collapse=", ")


df<-merge(df_antibody, df_lineage,by="X")
data<-df[,-(1:2)]
rownames(data)<-df$X

test_list <- c(colnames(data[, 1:(ncol(data))]))
test_list 

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$AML, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}


mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$Naïve_T, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}

colnames(mccf_final)

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$`0_HSPC`, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}


mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$AML, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}

colnames(mccf_final)

factor_list<-filter(mccf_final, normalized_mcc>0.5)
factor_list<-unique(factor_list$Factor)
mccf_plot<-filter(mccf_final, Factor%in%factor_list)


p<-ggplot(mccf_plot, aes(x = f1, y = normalized_mcc, fill=Factor))+
  geom_point(size = 2.5, shape = 21)+ylim(0,1)+xlim(0,1)+
  theme_minimal()+
  geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
  theme(plot.title = ggplot2::element_text(hjust = 0.5))+
  ylab("Normalized MCC")+xlab("F1 Score")+
  theme(axis.text.x =element_text(face="bold", 
                                  size=10),
        axis.text.y = element_text(face="bold", 
                                   size=10))+
  coord_equal(ratio = 1)

p

#############################################################
df_GE <- read.delim2("/oak/stanford/groups/cgawad/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/ALSF_RNA-seq_expr_AML_0.2.csv", 
                     sep = ",", header = TRUE)
df_GE[,2:(ncol(df_GE))]<-lapply(df_GE[,2:(ncol(df_GE))], as.numeric)
df_obs <- read.delim2("/oak/stanford/groups/cgawad/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/ALSF_RNA-seq_Combo_obs.csv",
                      sep = ",", header = TRUE)
colnames(df_GE)
df_GE <- df_GE %>% mutate(across(c(NOC2L:VBP1), ~ ifelse(.>0,1,0)))
#write.csv(df_GE, "ALSF_RNA-seq_expr_HSPC_Absent_binary.csv")

colnames(df_lineage)

df_lineage<-subset(df_obs, select = c("X","lineage"))
df_lineage<-df_lineage%>%mutate(Value=1)
df_lineage<-spread(df_lineage, key="lineage", value = "Value")
df_lineage <- df_lineage%>%mutate(across(c(`0_HSPC`:`T`), ~replace(., is.na(.), 0)))

df<-merge(df_GE, df_lineage,by="X")


data<-df[,-1]
rownames(data)<-df[,1]

test_list <- c(colnames(data[, 1:(ncol(data))]))
test_list 

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$AML, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}
write.csv(mccf_final, "ALSF_AML_MCCF_AML_AML_0.2.csv")
colnames(mccf_final)
factor_list<-filter(mccf_final, normalized_mcc>0.6)
factor_list<-filter(factor_list, f1>0.5)
factor_list<-unique(factor_list$Factor)

factor_list<-c("ENO1",
               "ITGA4",
               "CD74",
               "EGFL7",
               "CD44",
               "CD63",
               "SPN",
               "TYROBP",
               "MSN",
               'CSF3R')
factor_list<-c("TNFRSF14","CD34","ITGA6","TFPI", "CD38","PROM1","SLC38A1","ANPEP", "MBP","NECTIN2")

#"ARMH1, IGFBP7, SERPINB1, SOX4, MYB, CDK6, ETV6, GIHCG, PRSS57, CLEC11A, AML"
 
mccf_plot<-filter(mccf_final, Factor%in%factor_list)

p<-ggplot(mccf_plot, aes(x = f1, y = normalized_mcc, fill=Factor))+
  geom_point(size = 2.5, shape = 21)+ylim(0,1)+xlim(0,1)+
  theme_minimal()+
  scale_fill_manual(values=met.brewer("Signac", 13))+
  #geom_text(aes(label=ifelse(normalized_mcc>0.82,as.character(thresholds),'')),hjust=0,vjust=0)+
  geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
  theme(plot.title = ggplot2::element_text(hjust = 0.5))+
  ylab("Normalized MCC")+xlab("F1 Score")+
  theme(axis.text.x =element_text(face="bold", 
                                  size=12),
        axis.text.y = element_text(face="bold", 
                                   size=12))+
  coord_equal(ratio = 1)

p

paste(shQuote(factor_list), collapse=", ")

paste(factor_list, collapse=", ")

cellsurface<-read.table("GOCC_CELL_SURFACE.v2023.1.Hs.gmt", header=FALSE)
cellsurface_list<-as.list(cellsurface[1,c(3:931)])
intersect(df_2$term,cellsurface_list)

colnames(df_obs)

df_lineage<-subset(df_obs, select = c("X","Prognosis"))
df_lineage<-filter(df_lineage, !Prognosis=="0_HealthyBM")
df_lineage<-df_lineage%>%mutate(Prognosis_bi=1)%>%
  mutate(Prognosis_bi=ifelse(Prognosis=="Alive",0,Prognosis_bi))
colnames(df_lineage)

df<-merge(df_GE, df_lineage[,c(1,3)],by="X")

df$Prognosis_bi
data<-df[,-1]
rownames(data)<-df[,1]

test_list <- c(colnames(data[, 1:(ncol(data))]))
test_list 

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$Prognosis_bi , data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}
write.csv(mccf_final, "ALSF_AML_MCCF_AML_AML_0.2_prognosis.csv")
colnames(mccf_final)
factor_list<-filter(mccf_final, normalized_mcc>0.5)
factor_list<-filter(factor_list, f1>0.4)
factor_list<-unique(factor_list$Factor)
paste(shQuote(factor_list), collapse=", ")
###############################################################################
df_KMES <- read.delim2("ALSF_AML_ALL_cell_LSC_Module_MEs_pct_matrix.csv", sep = ",", header = TRUE)
colnames(df_KMES)
df_KMES<-df_KMES %>% separate(col="Sample", sep = "-", c("Sample", "Drop1"))
df_KMES[,1:(ncol(df_KMES)-2)]<-lapply(df_KMES[,1:(ncol(df_KMES)-2)], as.numeric)

df<-merge(df_clic,df_KMES[,-19],by='Sample')


data<-df[,-1]
rownames(data)<-df$Sample



test_list <- c(colnames(data[, 1:(ncol(data))]))
test_list 

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$Prognosis, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}

colnames(mccf_final)
factor_list<-filter(mccf_final, normalized_mcc>0.7)
factor_list<-unique(factor_list$Factor)
factor_list<-c(
  'MRD'   )

mccf_plot<-filter(mccf_final, Factor%in%factor_list)


p<-ggplot(mccf_plot, aes(x = f1, y = normalized_mcc, fill=Factor))+
  geom_point(size = 2.5, shape = 21)+ylim(0,1)+xlim(0,1)+
  theme_minimal()+
  scale_fill_manual(values=met.brewer("Signac", 10))+
  #geom_text(aes(label=ifelse(normalized_mcc>0.82,as.character(thresholds),'')),hjust=0,vjust=0)+
  geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
  theme(plot.title = ggplot2::element_text(hjust = 0.5))+
  ylab("Normalized MCC")+xlab("F1 Score")+
  theme(axis.text.x =element_text(face="bold", 
                                  size=12),
        axis.text.y = element_text(face="bold", 
                                   size=12))+
  coord_equal(ratio = 1)

p

autoplot(mccf1(data$Prognosis, data$GATA5_...))
###############################################################
df_obs <- read.delim2("ALSF_AML_Combo_obs.csv", sep = ",", header = TRUE)
colnames(df_obs)
df<-subset(df_obs, select = c("X","Prognosis", "CD19_pos",                    
                              "CD274_pos" ,                   "CD3_pos" ,                    
                              "CD33_pos"  ,                   "CD90_pos" ,                   
                              "CD10_pos"  ,                   "CD45RA_pos" ,                 
                              "CD123_pos"  ,                  "CD7_pos" ,                    
                              "CD49f_pos"  ,                  "CD25_pos" ,                   
                              "CD279_pos" ,                   "CD32_pos" ,                   
                              "CD152_pos"  ,                  "CD366_pos" ,                  
                              "CD235ab_pos",                  "CD127_pos" ,                  
                              "CD71_pos",                     "CD36_pos",                    
                              "CD133_pos"   ))

df[,3:(ncol(df))]<-lapply(df[,3:(ncol(df))], as.numeric)
df<-filter(df, !Prognosis=="0_HealthyBM")
df<-df%>%mutate(Prognosis_bi=1)%>%
  mutate(Prognosis_bi=ifelse(Prognosis=="Alive",0,Prognosis_bi))
colnames(df)

df<-merge(df,df_lineage,by="X")

data<-df[,-c(1,2)]
rownames(data)<-df[,1]

test_list <- c(colnames(data[, 1:(ncol(data))]))
test_list 

mccf_final<-data.frame(matrix(ncol = 4, nrow = 0))

for(t in test_list) {       # for-loop over columns
  mccf_df<-mccf1(data$AML, data[ , t])
  df <- as.data.frame(do.call(cbind, mccf_df))%>%mutate(Factor=paste0(t))
  mccf_final<-rbind(mccf_final,df)
}

colnames(mccf_final)
factor_list<-filter(mccf_final, normalized_mcc>0.5)
factor_list<-filter(factor_list, f1>0.5)
factor_list<-unique(factor_list$Factor)

paste(shQuote(factor_list), collapse=", ")

paste(factor_list, collapse=", ")

mccf_plot<-filter(mccf_final, Factor%in%factor_list)

p<-ggplot(mccf_plot, aes(x = f1, y = normalized_mcc, fill=Factor))+
  geom_point(size = 2.5, shape = 21)+ylim(0,1)+xlim(0,1)+
  theme_minimal()+
  scale_fill_manual(values=met.brewer("Signac", 8))+
  #geom_text(aes(label=ifelse(normalized_mcc>0.82,as.character(thresholds),'')),hjust=0,vjust=0)+
  geom_hline(yintercept=0.5, linetype="dashed", color = "red")+
  theme(plot.title = ggplot2::element_text(hjust = 0.5))+
  ylab("Normalized MCC")+xlab("F1 Score")+
  theme(axis.text.x =element_text(face="bold", 
                                  size=12),
        axis.text.y = element_text(face="bold", 
                                   size=12))+
  coord_equal(ratio = 1)

p
