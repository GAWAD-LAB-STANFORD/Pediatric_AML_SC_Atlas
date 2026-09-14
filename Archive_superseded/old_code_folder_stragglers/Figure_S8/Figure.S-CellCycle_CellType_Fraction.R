# ====================================================================
# Figure.S-CellCycle_CellType_Fraction.R  |  Figure S8 (composition + cell cycle)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : AML-sample-infor.csv, AML_scRNA_Cell_phase.csv, AML_scRNA_Cell_summary__Cell_Type_by_Sample.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript Figure.S-CellCycle_CellType_Fraction.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

rm(list=ls())
suppressPackageStartupMessages({
  library(tidyverse)
})
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggsci)
library(scales)
library(MetBrewer)
library(ggsignif)
library(grid)



setwd("~/Downloads/ALSF_AML_plot")

df_clic <- read.delim2("AML-sample-infor.csv", sep = ",", header = TRUE)
###########################################################
df <- read.delim2("AML_scRNA_Cell_phase.csv", sep = ",", header = TRUE)
df<-df%>% separate(col="Sample", sep = "-", c("Sample", "Drop1"))
df_L<-df[,-c(2,6)]
df_L<-gather(df_L, key = 'Cell Phase', value = 'Counts',-Sample)

library(RColorBrewer)
coul <- brewer.pal(3, "Pastel2") 

ggp <- ggplot(df_L,           # Create ggplot2 plot scaled to 1.00
              aes(x = Sample,
                  y = Counts,
                  fill = `Cell Phase`)) +
  geom_bar(position = "fill", stat = "identity")+
  scale_fill_brewer(palette = "Set2", direction = 1) +
  ylab('Percentage')+xlab("")+
  theme_minimal()+
  theme(axis.text.x = element_text(size=8, hjust=0.8,angle = 90),
        axis.text.y = element_text(size=10))
ggp  
#############################################
colnames(df)
df_P <- df%>%mutate(across(c(G1:sum), function(x) as.numeric(as.character(x))))
df_P_F<-df_P%>%mutate(across(c(G1:S), function(x) x/df_P$sum))
df_P_F<-df_P_F[,-c(2,6)]

df_P_F_plot<-merge(df_P_F, df_clic[,c(2,5)], all=T)
df_P_F_plot<-df_P_F_plot%>%mutate(Cytogenetic=ifelse(Sample%in%c("0_HealthyBM1","0_HealthyBM2"), 
                                                     'HBM', Cytogenetic))

df_P_F_plot<-gather(df_P_F_plot, key='Cell Phase', value='Percentage',-Cytogenetic,-Sample)

level_order=c("HBM","CN","del7q" ,"NUP98/NSD1","MLLr",
         "Tri(8)/MLLr","t(2;3)(p15;q26.2)",
         "t(7;14)(q21;q32)","BCR/ABL","Tri(15)" ,
         "PML/RARA", "CBFB/MYH11", "RUNX1/RUNX1T1",'MYB/GATA1','Tri(8)')

df_plot_p<- ggplot(df_P_F_plot,aes(x=factor(Cytogenetic,,levels = level_order),
                       y=Percentage,
                       color=Cytogenetic,
                      fill=Cytogenetic))+ 
  geom_boxplot() +
  theme_minimal()+
  ylab('Percentage')+xlab("")+
  scale_fill_manual(values=met.brewer("Egypt", 15))+
  theme(panel.grid = element_line(linetype = 15))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  expand_limits(y = c(0,1))+
  theme(axis.text.x = element_text(size=8, hjust=0.8,angle = 90),
        axis.text.y = element_text(size=10))+
  facet_grid(.~factor(`Cell Phase`),
             scales = "free_x")+
 theme(strip.text = element_text(size = 8))
#+
# geom_signif(comparisons = list(c("HBM","CN" ),
#                                c("HBM","MLLr"),
 #                                c("HBM","PML/RARA" ),
  #                               c("HBM","CBFB/MYH11"),
   #         map_signif_level = TRUE,
    #         step_increase = 0.1,
    #        test = 'wilcox.test',
     #         test.args=list(alternative = "two.sided", var.equal = FALSE, paired=FALSE),
      #        y_position = 1,tip_length = 0, vjust = 0.2)


print(df_plot_p)
# Finished line plot
  theme(axis.text.x =element_blank(),
        axis.text.y = element_text(face="bold", 
                                   size=10))

dev.off()

df_P_F_plot<-merge(df_P_F, df_clic[,c(2,35)], all=T)
df_P_F_plot<-df_P_F_plot%>%mutate(Alive.deceased=ifelse(Sample%in%c("0_HealthyBM1","0_HealthyBM2"), 
                                                     'HBM',Alive.deceased))
df_P_F_plot<-filter(df_P_F_plot, !Sample=='AML4363')
df_P_F_plot<-gather(df_P_F_plot, key='Cell Phase', value='Percentage',-Alive.deceased,-Sample)

df_plot_p<- ggplot(df_P_F_plot,aes(x=Alive.deceased,
                                   y=Percentage,
                                   color=Alive.deceased,
                                   fill=Alive.deceased))+ 
  geom_boxplot() +
  theme_minimal()+
  ylab('Percentage')+xlab("")+
  scale_fill_manual(values=met.brewer("Egypt", 3))+
  theme(panel.grid = element_line(linetype = 4))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  expand_limits(y = c(0,1))+
  theme(axis.text.x = element_text(size=8, hjust=0.8,angle = 90),
        axis.text.y = element_text(size=10))+
  facet_grid(.~factor(`Cell Phase`),
             scales = "free_x")+
  theme(strip.text = element_text(size = 8))+
  geom_signif(comparisons = list(c("HBM","Alive"),
                              c("HBM","Deceased"),
                              c("Alive","Deceased")),
              map_signif_level = TRUE,
              step_increase = 0.06,
              test = 't.test',
              size = 0.5,
              textsize = 2,
              test.args=list(alternative = "two.sided", var.equal = FALSE, paired=FALSE),
              y_position = 1,tip_length = 0, vjust = 0.2)


print(df_plot_p)

dev.off()


#############################
df_cell <- read.delim2("Figure_S12/AML_scRNA_Cell_summary__Cell_Type_by_Sample.csv", sep = ",", header = TRUE)
colnames(df_cell)
df_cell<-df_cell%>% separate(col="Sample", sep = "-", c("Sample", "Drop1"))
df_cell_L<-df_cell[,-c(2,ncol(df_cell))]
df_cell_L<-gather(df_cell_L, key = 'Cell Type', value = 'Counts',-Sample)


#Scanpy color for cell types
col_list=c('0_HSPC'= '#023fa5',
           'AML'='#7d87b9',
           'AML.B'= '#bec1d4',
           'AML.CD14'= '#d6bcc0',
           'AML.CD1C'= '#bb7784',
           'AML.CD4T'= '#8e063b',
           'AML.CTL'= '#4a6fe3',
           'AML.Ery'= '#8595e1',
           'AML.MKI67'= '#b5bbe3',
           'AML.NK'= '#e6afb9',
           'AML.Naïve_CD8T'= '#e07b91',
           'AML.PCNA'= '#d33f6a',
           'Activated_CD4T'= '#11c638',
           'CD14_Monocyte'= '#8dd593',
           'CD16_Monocyte'= '#c6dec7',
           'CD20.B'= '#ead3c6',
           'CD34.ProB'= '#f0b98d',
           'CTL'= '#ef9708',
           'Erythrocytes'= '#0fcfc0',
           'Macrophage'= '#9cded6',
           'Myeloid_Pro'= '#d5eae7',
           'NK'= '#f3e1eb',
           'Naïve_CD4T'= '#f6c4e1',
           'Naïve_CD8T'= '#f79cd4',
           'PlasmaB'= '#7f7f7f',
           'PreB'= '#c7c7c7',
           'ProB'= '#1ce6ff',
           'mDC'= '#336600',
           'pDC'= '#99adc0')

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


df_Cell_P<- df_cell%>%mutate(across(c(X0_HSPC:Cell_sum), function(x) as.numeric(as.character(x))))
df_Cell_P<-df_Cell_P%>%mutate(across(c(X0_HSPC:pDC), function(x) x/df_Cell_P$Cell_sum))
df_Cell_P<-df_Cell_P%>%mutate(mature_B_sum=0,CTL_sum=0,
                              NK_sum=0,CD4_sum=0,
                              Naïve_CD8_sum=0,
                              CD14_Mono_sum=0,
                              Ery_sum=0,
                              mDC_sum=0)
df_Cell_P<-df_Cell_P%>%mutate(mature_B_sum=df_Cell_P$AML.B+df_Cell_P$CD20.B)
df_Cell_P<-df_Cell_P%>%mutate(CTL_sum=df_Cell_P$AML.CTL+df_Cell_P$CTL)
df_Cell_P<-df_Cell_P%>%mutate(NK_sum=df_Cell_P$AML.NK+df_Cell_P$NK)
df_Cell_P<-df_Cell_P%>%mutate(CD4_sum=df_Cell_P$AML.CD4T+df_Cell_P$Naïve_CD4T+df_Cell_P$Activated_CD4T)
df_Cell_P<-df_Cell_P%>%mutate(Naïve_CD8_sum=df_Cell_P$AML.Naïve_CD8T+df_Cell_P$Naïve_CD8T)
df_Cell_P<-df_Cell_P%>%mutate(CD14_Mono_sum=df_Cell_P$AML.CD14+df_Cell_P$CD14_Monocyte)
df_Cell_P<-df_Cell_P%>%mutate(Ery_sum=df_Cell_P$AML.Ery+df_Cell_P$Erythrocytes)
df_Cell_P<-df_Cell_P%>%mutate(mDC_sum=df_Cell_P$AML.CD1C+df_Cell_P$mDC)



df_Cell_lineage<-df_Cell_P
colnames(df_Cell_lineage)


df_Corela2 <- df_Cell_lineage[,3:31]
rownames(df_Corela2) <- df_Cell_lineage[,1]

Corela_df<-round(cor(df_Corela2),digits = 2 )

library(corrplot)
testRes<-cor.mtest(df_Corela2,conf.level = 0.95)
M = cor(df_Corela2)

col=colorRampPalette(met.brewer("Signac", 2))(200)

corrplot(M, p.mat = testRes$p, method = 'circle',insig='blank',
         type ='lower',
         diag = FALSE,
         tl.cex = 0.8,addCoef.col = 1, tl.col = "brown",   # Change font size of correlation coefficients
         number.cex = 0.6,
         tl.srt = 45,
         order = 'AOE',
)$corrPos -> p1



df_B_lineage<-subset(df_Cell_lineage,select=c("Sample","CD34.ProB" ,
                                                 "ProB" ,"PreB", "PlasmaB",
                                              "CD20.B", "AML.B"))

df_plot<-df_Cell_lineage
df_plot$`Sample Type`<-'AML'

df_plot<-df_plot%>%
  mutate(`Sample Type`=ifelse(Sample%in%c("0_HealthyBM1","0_HealthyBM2"), 
                                                          'HBM',`Sample Type`))

ggplot(df_plot) +
  geom_line(aes(x = Sample, y=CD34.ProB,color="CD34.ProB"),group = 1) +
  geom_line(aes(x = Sample, y=ProB, color= 'ProB', group = 1)) +
  geom_line(aes(x = Sample, y=PreB, color= 'PreB',group = 1)) +
  geom_line(aes(x = Sample, y=CD20.B, color= 'CD20.B',group = 1)) +
  geom_line(aes(x = Sample, y=PlasmaB, color= 'PlasmaB',group = 1)) +
  geom_line(aes(x = Sample, y=AML.B, color= 'AML.B',group = 1)) +
  labs(x = "Sample",
       y = "Percentage",
       color = col_list) +
  scale_color_manual(values = col_list)+
  theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



df_plot<-df_plot[,-c(2,32)]%>%gather(key='Cell_Type', value='Percentage', -Sample,-`Sample Type`)


df_plot_set<-df_plot[-grep('AML',df_plot$Cell_Type),]


unique(df_plot_set$Cell_Type)

cell_type_level<-c("X0_HSPC"  ,    "Myeloid_Pro",   "Erythrocytes" ,
                   "CD14_Monocyte",  "CD16_Monocyte" ,"Macrophage",
                   "mDC"   ,         "pDC" ,   
                   "CD34.ProB" ,         "ProB"  ,    "PreB" ,  "CD20.B"  ,   "PlasmaB" ,             
                   "Activated_CD4T","Naïve_CD4T" , "Naïve_CD8T" ,"CTL"  ,    
                   "NK" , 
                   'CD14_Mono_sum',
                   'Ery_sum',
                   'mDC_sum',
                   "mature_B_sum", 
                   "CD4_sum","Naïve_CD8_sum","CTL_sum","NK_sum")

df_plot_p<- ggplot(df_plot_set,aes(x=factor(`Sample Type`,
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
  facet_wrap(.~factor(`Cell_Type`,
                      levels = cell_type_level),ncol = 4,
             scales = "free_y",
             #strip.position="right"
             )+
  theme(strip.text = element_text(size = 8),
        panel.spacing = unit(1, "lines"))+
  geom_signif(comparisons = list(c("HBM","AML")),
              map_signif_level = TRUE,
              test = 'wilcox.test',
              size = 0.3,
              textsize = 4,
              #test.args=list(alternative = "two.sided", var.equal = FALSE, paired=FALSE),
              y_position =0.005,
              tip_length = 0, 
              vjust = 0.2)



print(df_plot_p)

dev.off()



colnames(df_Cell_lineage)

df_T_lineage<-subset(df_Cell_lineage,select=c("Sample","AML.CTL",
                                                 "CTL" ,"AML.CD4T",
                                                 "Naïve_CD4T",
                                                 "Naïve_CD8T",
                                                 'AML.Naïve_CD8T'))

ggplot(df_T_lineage) +
  geom_line(aes(x = Sample, y=AML.CTL,color="AML.CTL"),group = 1) +
  geom_line(aes(x = Sample, y=CTL, color= 'CTL', group = 1)) +
  geom_line(aes(x = Sample, y=AML.CD4T,color="AML.CD4T"),group = 1) +
  geom_line(aes(x = Sample, y=Naïve_CD4T, color= "Naïve_CD4T", group = 1)) +
  geom_line(aes(x = Sample, y=AML.Naïve_CD8T,color="AML.Naïve_CD8T"),group = 1) +
  geom_line(aes(x = Sample, y=Naïve_CD8T, color= "Naïve_CD8T", group = 1)) +
  labs(x = "Sample",
       y = "Percentage",
       color = col_list) +
  scale_color_manual(values = colors)+
  theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


df_leukemia<-filter(df_Cell_lineage[,-c(2,ncol(df_Cell_lineage))],
                    Sample%in%c("0_HealthyBM1","0_HealthyBM2",
                               "AML4035",'AML4232','AML4062','AML4127'))
df_leukemia_2<-t(df_leukemia)

colnames(df_leukemia_2)<-df_leukemia[,1]
data<-df_leukemia_2[-1,]
data<-as.data.frame(data)
data$`Cell Type`<-rownames(data)
colnames(data)
data <- data%>%mutate(across(c(`0_HealthyBM1`:`AML4232`), 
                             function(x) as.numeric(as.character(x))))

col_list=c('0_HealthyBM1'='#98df8a',
           '0_HealthyBM2'= '#c5b0d5',
           'AML4035' = '#ff9896',
           'AML4232'='#aa40fc',
           'AML4062'= '#aec7e8',
           'AML4127'= '#17becf')

ggplot(data) +
  geom_line(aes(x = `Cell Type`, y=`0_HealthyBM1`,color= '0_HealthyBM1'),group = 1)+
  geom_line(aes(x = `Cell Type`, y=`0_HealthyBM2`, color= '0_HealthyBM2'), group = 1) +
  geom_line(aes(x =`Cell Type`, y=AML4232, color= "AML4232"), group = 1) +
  geom_line(aes(x = `Cell Type`, y=AML4062,color="AML4062"),group = 1)+
  labs(x = "Cell Type",
       y = "Percentage",
       color=col_list) +
  scale_fill_manual(values = col_list)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

##########################

df_cell <- read.delim2("Figure_S12/AML_scRNA_Cell_summary__Cell_Type_by_Sample.csv", sep = ",", header = TRUE)
colnames(df_cell)
df_cell<-df_cell%>% separate(col="Sample", sep = "-", c("Sample", "Drop1"))


df_Cell_P<- df_cell%>%mutate(across(c(X0_HSPC:Cell_sum), function(x) as.numeric(as.character(x))))

colnames(df_Cell_P)

column_list<-c( "AML","AML.MKI67","AML.PCNA",'Cell_sum')

df_Cell_P<-df_Cell_P[, -which(names(df_Cell_P) %in%column_list)]

df_Cell_P$Cell_sum <- apply(df_Cell_P[,3:28], 1, sum)
df_Cell_P<-df_Cell_P%>%mutate(across(c(X0_HSPC:pDC), function(x) x/df_Cell_P$Cell_sum))
df_Cell_P<-df_Cell_P%>%mutate(mature_B_sum=0,CTL_sum=0,
                              NK_sum=0,CD4_sum=0,
                              Naïve_CD8_sum=0,
                              CD14_Mono_sum=0,
                              Ery_sum=0,
                              mDC_sum=0)
df_Cell_P<-df_Cell_P%>%mutate(mature_B_sum=df_Cell_P$AML.B+df_Cell_P$CD20.B)
df_Cell_P<-df_Cell_P%>%mutate(CTL_sum=df_Cell_P$AML.CTL+df_Cell_P$CTL)
df_Cell_P<-df_Cell_P%>%mutate(NK_sum=df_Cell_P$AML.NK+df_Cell_P$NK)
df_Cell_P<-df_Cell_P%>%mutate(CD4_sum=df_Cell_P$AML.CD4T+df_Cell_P$Naïve_CD4T+df_Cell_P$Activated_CD4T)
df_Cell_P<-df_Cell_P%>%mutate(Naïve_CD8_sum=df_Cell_P$AML.Naïve_CD8T+df_Cell_P$Naïve_CD8T)
df_Cell_P<-df_Cell_P%>%mutate(CD14_Mono_sum=df_Cell_P$AML.CD14+df_Cell_P$CD14_Monocyte)
df_Cell_P<-df_Cell_P%>%mutate(Ery_sum=df_Cell_P$AML.Ery+df_Cell_P$Erythrocytes)
df_Cell_P<-df_Cell_P%>%mutate(mDC_sum=df_Cell_P$AML.CD1C+df_Cell_P$mDC)


df_plot<-df_Cell_P
df_plot$`Sample Type`<-'AML'

df_plot<-df_plot%>%
  mutate(`Sample Type`=ifelse(Sample%in%c("0_HealthyBM1","0_HealthyBM2"), 
                              'HBM',`Sample Type`))

df_plot<-df_plot[,-c(2,29)]%>%gather(key='Cell_Type', value='Percentage', -Sample,-`Sample Type`)

unique(df_plot$Cell_Type)

cell_type_level<-c("X0_HSPC"  ,    "Myeloid_Pro",   "Erythrocytes" ,
                   "CD14_Monocyte",  "CD16_Monocyte" ,"Macrophage",
                   "mDC"   ,         "pDC" ,   
                   "CD34.ProB" ,         "ProB"  ,    "PreB" ,  "CD20.B"  ,   "PlasmaB" ,             
                   "Activated_CD4T","Naïve_CD4T" , "Naïve_CD8T" ,"CTL"  ,    
                   "NK" , 
                   'CD14_Mono_sum',
                   'Ery_sum',
                   'mDC_sum',
                   "mature_B_sum", 
                   "CD4_sum","Naïve_CD8_sum","CTL_sum","NK_sum")
df_plot_set<-df_plot[df_plot$Cell_Type%in%cell_type_level,]
df_plot_p<- ggplot(df_plot_set,aes(x=factor(`Sample Type`,
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
  facet_wrap(.~factor(`Cell_Type`,
                      levels = cell_type_level),ncol = 4,
             scales = "free_y",
             #strip.position="right"
  )+
  theme(strip.text = element_text(size = 8),
        panel.spacing = unit(1, "lines"))+
  geom_signif(comparisons = list(c("HBM","AML")),
              map_signif_level = TRUE,
              test = 'wilcox.test',
              size = 0.3,
              textsize = 4,
              #test.args=list(alternative = "two.sided", var.equal = FALSE, paired=FALSE),
              y_position =0.005,
              tip_length = 0, 
              vjust = 0.2)



print(df_plot_p)

dev.off()

