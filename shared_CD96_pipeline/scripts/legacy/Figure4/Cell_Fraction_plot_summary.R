# ====================================================================
# Cell_Fraction_plot_summary.R  |  Figure 4 (T/NK states + survival)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : AML_scRNA_T_Cell_Cell_Type.csv, AML-sample-infor.csv, Combo_gene_fraction_bySample.csv, AML_scRNA_Cell_summary__PAC_anno_by_Sample.csv, Gene_feature_score_Combo_enriched.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript Cell_Fraction_plot_summary.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
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


setwd("~/Downloads/ALSF_AML_plot")
##################################################
df_clic <- read.delim2("AML-sample-infor.csv", sep = ",", header = TRUE)
df_clic <-df_clic %>%mutate(Prognosis=0)%>%mutate(Prognosis=ifelse(Alive.deceased=='Alive',0,Prognosis))%>%
  mutate(Prognosis=ifelse(Alive.deceased=='Deceased',1,Prognosis))

df_clic$MRD<-df_clic$MRD...
df_clic$Blast<-df_clic$Blast..
df_clic<-subset(df_clic,select=c("Sample","FAB_2","Cytogenetic","Alive","Relapsed","Remission","Prognosis","MRD","Blast" ))
colnames(df_clic)
df_clic <- df_clic%>%mutate(across(c(Alive:Blast), function(x) as.numeric(as.character(x))))
df_clic <- df_clic%>%mutate(across(c(Alive:Blast), ~replace(., is.na(.), 0)))
df_clic <- filter(df_clic, !Sample=="AML4363")
###################################################
df_cell <- read.delim2("ALSF AML Figure5/AML_scRNA_T_Cell_Cell_Type.csv", sep = ",", header = TRUE)
colnames(df_cell)
df_cell$Sample

df_cell<-df_cell%>% separate(col="Sample", sep = "-", c("Sample", "Drop1"))
df_cell_L<-df_cell[,-c(2,ncol(df_cell))]
df_cell_L<-gather(df_cell_L, key = 'Cell Type', value = 'Counts',-Sample)
unique(df_cell_L$`Cell Type`)

#Scanpy color for cell types
col_list=c('AML.Activated.CD4T'='#1f77b4',
           'AML.GZMB.CD8T'='#ff7f0e',
           'AML.GZMK.CD8T'='#279e68',
           'AML.MAIT'='#d62728',
           'AML.NK'='#aa40fc',
           'AML.Naïve.CD4T'='#8c564b',
           'AML.Naïve.CD8T'='#e377c2',
           'Activated.CD4T'='#b5bd61',
           'Effector.memory.CD8T'='#17becf',
           'GZMB.CD8T'='#aec7e8',
           'GZMB.DNT'='#ffbb78',
           'GZMB.NK'='#98df8a',
           'GZMK.CD8T'='#ff9896',
           'GZMK.NK'='#c5b0d5',
           'MAIT'='#c49c94',
           'Naïve.CD4T'='#f7b6d2',
           'Naïve.CD8T'='#dbdb8d'
          )

ggp <- ggplot(df_cell_L,           # Create ggplot2 plot scaled to 1.00
              aes(x = Sample,
                  y = Counts,
                  fill = `Cell Type`)) +
  geom_bar(position = "fill", stat = "identity")+
  ylab('Percentage')+xlab("")+
  scale_fill_manual(values = col_list)+
  theme_minimal()+
  theme(axis.text.x = element_text(size=8, hjust=0.8,angle = 90),
        axis.text.y = element_text(size=10))
ggp  


df_Cell_P<- df_cell%>%mutate(across(c(AML.Activated.CD4T:sum), function(x) as.numeric(as.character(x))))
df_Cell_P<-df_Cell_P%>%mutate(across(c(AML.Activated.CD4T:Naïve.CD8T), function(x) x/df_Cell_P$sum))
df_Cell_P<-df_Cell_P%>%mutate(NaïveCD4_sum=0,
                              NaïveCD8_sum=0,
                              Activated.CD4_sum=0,
                              GZMBCD8_sum=0,
                              GZMKCD8_sum=0,
                              MAIT_sum=0,
                              GZMBNK_sum=0)
df_Cell_P<-df_Cell_P%>%mutate(NaïveCD4_sum=df_Cell_P$AML.Naïve.CD4T+df_Cell_P$Naïve.CD4T)
df_Cell_P<-df_Cell_P%>%mutate(NaïveCD8_sum=df_Cell_P$AML.Naïve.CD8T+df_Cell_P$Naïve.CD8T)
df_Cell_P<-df_Cell_P%>%mutate(Activated.CD4_sum=df_Cell_P$Activated.CD4T+df_Cell_P$AML.Activated.CD4T)
df_Cell_P<-df_Cell_P%>%mutate(GZMBCD8_sum=df_Cell_P$AML.GZMB.CD8T+df_Cell_P$GZMB.CD8T)
df_Cell_P<-df_Cell_P%>%mutate(GZMKCD8_sum=df_Cell_P$AML.GZMK.CD8T+df_Cell_P$GZMK.CD8T)
df_Cell_P<-df_Cell_P%>%mutate(MAIT_sum=df_Cell_P$AML.MAIT+df_Cell_P$MAIT)
df_Cell_P<-df_Cell_P%>%mutate(GZMBNK_sum=df_Cell_P$AML.NK+df_Cell_P$GZMB.NK)
df_Cell_P<-df_Cell_P%>%mutate(GZMKNK_sum=df_Cell_P$GZMK.NK)
df_Cell_lineage<-df_Cell_P
colnames(df_Cell_lineage)




df_plot<-df_Cell_lineage
df_plot$`Sample Type`<-'AML'

df_plot<-df_plot%>%
  mutate(`Sample Type`=ifelse(Sample%in%c("0_HealthyBM1","0_HealthyBM2"), 
                              'HBM',`Sample Type`))


df_plot<-df_plot[,-c(2,20)]%>%gather(key='Cell_Type', value='Percentage', -Sample,-`Sample Type`)


df_plot_set<-df_plot[-grep('AML',df_plot$Cell_Type),]


unique(df_plot_set$Cell_Type)

df_plot_p<- ggplot(df_plot,aes(x=factor(`Sample Type`,
                                            levels = c('HBM','AML')),
                                   y=Percentage,
                                   color=`Sample Type`,
                                   fill=`Sample Type`))+ 
  geom_boxplot(width=0.3) +
  theme_minimal()+
  ylab('Percentage')+xlab("")+
  scale_fill_manual(values=met.brewer("Egypt", 3))+
  theme(panel.grid = element_line(linetype = 4))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  #expand_limits(y = c(0,0.05))+
  theme(axis.text.x = element_text(size=8, hjust=0.8,angle = 90),
        axis.text.y = element_text(size=10))+
  facet_wrap(.~factor(`Cell_Type`),ncol = 4,
             scales = "free_y",
             #strip.position="right"
  )+
  theme(strip.text = element_text(size = 8),
        panel.spacing = unit(1, "lines"))+
  geom_signif(comparisons = list(c("HBM","AML")),
              map_signif_level = TRUE,
              test = 't.test',
              size = 0.3,
              textsize = 4,
              #test.args=list(alternative = "two.sided", var.equal = FALSE, paired=FALSE),
              y_position =0.005,
              tip_length = 0, 
              vjust = 0.2)



print(df_plot_p)

dev.off()

####################################
df_PAC <- read.delim2("Figure_T_Cells/AML_scRNA_Cell_summary__PAC_anno_by_Sample.csv", sep = ",", header = TRUE)
  
colnames(df_PAC)


df_PAC<-df_PAC%>% separate(col="Sample", sep = "-", c("Sample", "Drop1"))
df_PAC_L<-df_PAC[,-c(2,ncol(df_PAC))]
df_PAC_L<-gather(df_PAC_L, key = 'Cell Type', value = 'Counts',-Sample)

df_PAC_L<-filter(df_PAC_L,Counts>10)

df_PAC_P<- df_PAC%>%mutate(across(c(X0_HSPC:Cell_sum), function(x) as.numeric(as.character(x))))
df_PAC_P<-df_PAC_P%>%mutate(across(c(X0_HSPC:pDC), function(x) x/df_Cell_P$sum))
df_PAC_P<-df_PAC_P%>% separate(col="Sample", sep = "-", c("Sample", "Drop1"))

df<-merge(df_PAC_P[,c(1,22:25,31:35)],df_Cell_lineage[,-2], by='Sample')

data<-df[,-c(1,28)]


data[,1:ncol(data)]<-lapply(data[,1:ncol(data)], as.numeric)
rownames(data) <- df[,1]

library(corrplot)
testRes<-cor.mtest(data,conf.level = 0.95)
M = cor(data)

col=colorRampPalette(met.brewer("Signac", 2))(200)

corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         diag = FALSE,
         tl.cex = 0.8,addCoef.col = 1,    # Change font size of correlation coefficients
         number.cex = 0.6,
         tl.srt = 45,
         order = 'alpha',
         col=colorRampPalette(met.brewer("Archambault", 10))(10))$corrPos -> p1

corrplot(M, p.mat = testRes$p, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.4, order = 'AOE', diag=FALSE)


corrplot(M, order = 'AOE', col = COL2('RdBu', 10))


text(p1$x, p1$y, round(p1$corr, 2))

corrplot2(
  data =df_Corela2,
  method = "pearson",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 25
)

corrplot(Corela_df,, method = 'circle', type = 'lower', insig='blank',
         addCoef.col ='black', number.cex = 0.8, order = 'AOE', diag=FALSE)

corrplot.mixed(Corela_df,lower = 'shade', upper = 'pie', order = 'hclust')

Corela_df<-round(cor(df_Corela2),digits = 2 )

######################################################################################################

df_GF <- read.delim2("Gene_feature_score_Combo_enriched.csv", sep = ",", header = TRUE)
colnames(df_GF)
df_GF<- df_GF%>%mutate(across(c(HSC_Score:LSC_EPPERT_Score), function(x) as.numeric(as.character(x))))
df_GF_anno<-unique(select(df_GF,c("SampleID","Prognosis")))
str(df_GF)

df_GF_mean <- df_GF %>% 
   group_by(SampleID) %>% 
   summarise(across(ends_with("_Score"), list(mean = mean, sd = sd), na.rm = TRUE, .names = "{col}_{fn}"))

df<-merge(df_GF_anno, df_GF_mean, by='SampleID')

df<-filter(df, !SampleID==c('AML4363',"0_HealthyBM1" ,"0_HealthyBM2"))

colnames(df)
level_order<-c('Alive','Deceased')

df_p<-ggplot(df_GF, aes(x=Prognosis,
                                 y=LSC_EPPERT_Score ,
                                 fill=Prognosis) )+
  geom_boxplot()+ 
  theme_minimal()+
  scale_fill_manual(values=met.brewer("Egypt", 3))+
  theme(panel.grid = element_line(linetype = 2))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  theme(axis.text.x = element_text( vjust = 0.5,hjust=0.6))+
  geom_signif(comparisons = list(c("Alive", "Deceased")),
              y_position = 0.4,vjust = 0.2)
print(df_p)
# Finished line plot
pdf(paste0("ALSF_AML_HBM_VS_AML_",t,".pdf"), width = 4, height =4)
#############################################################################
df_cell<-gather(df_cell, key = 'Cell_Type', value = 'Counts',-batchName)

ggp <- ggplot(df_cell,           # Create ggplot2 plot scaled to 1.00
              aes(x = batchName,
                  y = Counts,
                  fill = Cell_Type)) +
  geom_bar(position = "fill", stat = "identity")+
  theme_minimal()+
  theme(axis.text.x = element_text( size=8, hjust=0.8,angle = 90),
        axis.text.y = element_text(face="bold",size=12))
ggp 
###########################################################################################
df_clic <- read.delim2("AML-sample-infor.csv", sep = ",", header = TRUE)
df_clic <-df_clic %>%mutate(Prognosis=0)%>%mutate(Prognosis=ifelse(Alive.deceased=='Alive',1,Prognosis))%>%
  mutate(Prognosis=ifelse(Alive.deceased=='Deceased',0,Prognosis))
df_clic<-subset(df_clic,select=c(2,3,4,5,11))
df_gene<-read.delim2("Combo_gene_fraction_bySample.csv", sep = ",", header = TRUE)

genelist<-c('Sample','NPR3', 'ANGPT1', 'MECOM','ATP8B4','SELENOP','CD34','HOPX',
      'BAALC','HPGDS','SLC45A3','RHEX','CLIP3','TIMP3','DEPTOR','ZNF385D',
      'SPNS2','HOXA9','MYCN','BAALC','SPINK2','PROM1'
      #'CCNA1','AC125603.2', 'AC244502.1','MSLN','PITX1','GTSF1','RUNX1T1','ST18','GGT5','EGF','HOXB6','MYO18B','NDST3','KRT17','AGXT','IL9R','GJA1',
      #'CPA3','HGF','MPO','AZU1','TSPOAP1','TRH','CEBPE','PRTN3','ELANE','PLPPR3','CTSG','IGFBP2',
      #'CD96','CD7','LAT','CCL3','GNLY', 'PRF1', 'CD3E', 'BCL11B','TRBC1',
      #'CD9','TPSAB1'
      )
df_filtered<-subset(df_gene, select=genelist)
HBM<-c('HealthyBM1','HealthyBM2')


df_filtered<-filter(df_filtered,!df_filtered$Sample%in%HBM)

df_Corela<-merge(df_clic,df_cell_P[,-2],by = 'Sample')

df_Corela[,4:27]<-lapply(df_Corela[,4:27], as.numeric)

df_Corela2 <- df_Corela[,4:27]
rownames(df_Corela2) <- df_Corela[,1]

Corela_df<-round(cor(df_Corela2),digits = 2 )

corrplot2 <- function(data,
                      method = "pearson",
                      sig.level = 0.05,
                      order = "original",
                      diag = FALSE,
                      type = "upper",
                      tl.srt = 90,
                      number.font = 1,
                      number.cex = 1,
                      mar = c(0, 0, 0, 0)) {
library(corrplot)
  data_incomplete <- data
  data <- data[complete.cases(data), ]
  mat <- cor(data, method = method)
  cor.mtest <- function(mat, method) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], method = method)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  p.mat <- cor.mtest(data, method = method)
  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  corrplot(mat,
           method = "color", col = col(200), number.font = number.font,
           mar = mar, number.cex = number.cex,
           type = type, order = order,
           addCoef.col = "black", # add correlation coefficient
           tl.col = "black", tl.srt = tl.srt, # rotation of text labels
           # combine with significance level
           p.mat = p.mat, sig.level = sig.level, insig = "blank",
           # hide correlation coefficients on the diagonal
           diag = diag
  )
}

# edit from here
corrplot2(
  data = df_Corela2,
  method = "pearson",
  sig.level = 0.05,
  order = "original",
  diag = FALSE,
  type = "upper",
  tl.srt = 75
)

