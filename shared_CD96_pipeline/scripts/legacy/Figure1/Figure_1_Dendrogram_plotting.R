# ====================================================================
# Figure_1_Dendrogram_plotting.R  |  Figure 1E (per-sample dendrogram)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : AML-sample-infor.csv, AML_scRNA_Cell_summary_Cell_Type.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript Figure_1_Dendrogram_plotting.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================


rm(list=ls())
suppressPackageStartupMessages({
  library(tidyverse)
})
library(ggplot2)
library(reshape2)
library(dplyr)
require(flexmix)          # For model bases clustering
require(ggplot2)          # For plotting
require(FactoMineR)       # For performing MCA
require(vegan)            # For calculating jaccard distance
require(gplots)           # For plotting heatmaps
# For plotting heatmaps with annotations
require(reshape)          # For melting dataframe
require(adegenet)         # For performing the Directed Minimum Spanning tree for clones

require(igraph)           # For plotting of trees 
require(ggplot2)
require(grid)
require(scales)
require(reshape2)
require(gplots)
require(RColorBrewer)
#for clustering
require(vegan)
require(cluster)

#For tree Drawing
require(phangorn)
require(ape)
require(adegenet)
library(dplyr)
library(tibble)
library(dendextend)
library(ComplexHeatmap)
library(circlize)

detach(package:plyr)

setwd("~/Downloads/ALSF_AML_plot")
df_pca <- read.table("pca_Combo_AML_with_normal_by_Sample_for_dendrogram_repca.csv",sep = ",", header = TRUE)
df_pca[2:ncol(df_pca)]<- lapply(df_pca[2:ncol(df_pca)], FUN = function(y){as.numeric(y)})
df_pca<-df_pca%>% separate(col=1, sep = "-", c("Sample", "Drop1"))
df_pca[3:ncol(df_pca)]<- lapply(df_pca[3:ncol(df_pca)], FUN = function(y){as.numeric(y)})

###################################
df_cell <- read.delim2("AML_scRNA_Cell_summary_Cell_Type.csv", sep = ",", header = TRUE)
df_Cell_corr<-df_cell %>% separate(col="Sample", sep = "-", c("Sample", "Drop1"))

df_Cell_corr<-subset(df_Cell_corr,select=c(1,4,8,9,15,18,36))
df_Cell_corr <- df_Cell_corr%>%mutate(across(c(AML:AML.PCNA), function(x) as.numeric(as.character(x))))
df_Cell_corr$AML_sum <- apply(df_Cell_corr[,2:6], 1, sum)
df_Cell_corr$LSC_sum <- apply(df_Cell_corr[,3:4], 1, sum)
df_Cell_corr$CC_sum <- apply(df_Cell_corr[,5:6], 1, sum)
colnames(df_Cell_corr)
df_Cell_Blast<- subset(df_Cell_corr,select=c(1))
df_Cell_Blast$AML<- df_Cell_corr$AML_sum/df_Cell_corr$Cell_sum
df_Cell_Blast$AML_CD34<- df_Cell_corr$LSC_sum/df_Cell_corr$AML_sum
df_Cell_Blast$AML_CC<- df_Cell_corr$CC_sum/df_Cell_corr$AML_sum

percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...))
}

df_Cell_Blast<-df_Cell_Blast%>%mutate(across(c(AML:AML_CC), function(x) percent(x)))

############################################################################################
df_clic <- read.delim2("AML-sample-infor.csv", sep = ",", header = TRUE)
df_clic <-df_clic %>%mutate(Prognosis=0)%>%mutate(Prognosis=ifelse(Alive.deceased=='Alive',0,Prognosis))%>%
  mutate(Prognosis=ifelse(Alive.deceased=='Deceased',1,Prognosis))
df_clic$MRD<-df_clic$MRD...
df_clic$Blast<-df_clic$Blast..
#df_clic <-merge(df_clic, df_pca, by="Sample", all = T)
##AML4271 and AML882 don't have clear diagnosis about relapse.
df_clic<-df_clic%>%mutate(Relapsed=ifelse(Sample=="AML4271"| Sample=="AML882", "none" ,Relapsed))
colnames(df_clic)

df_clic_subset<-df_clic[,c('Sample','Cytogenetic')]
###################################################################
####################################################################
df<-merge(df_pca, df_clic, by='Sample', all=F)
mydata<-df[3:52]
rownames(mydata)<-df$Sample

#mydata<-df_pca[3:32]
#rownames(mydata)<-df_pca$Sample
rows.cor <- cor(t(mydata), use = "complete.obs", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor))
hclust.row$labels
########################################
df_cytoge<-df%>%dplyr::select(c("Sample", "Cytogenetic"))
df_cytoge<-df_cytoge%>%mutate(value='TRUE')
df_cytoge_w<-spread(df_cytoge, key='Cytogenetic', value='value',fill='none')
colnames(df)
df_mutate<-cbind(df_cytoge_w[1:15],df[57:81])
df_mutate[is.na(df_mutate)] <- 'none'
df_mutate_L<-gather(df_mutate,key='Genetic Variant', value='value',-Sample)

cols <- c('TRUE' = '#F4A582', 'none' = 'white')
ggplot(df_mutate_L, aes(x= Sample, y = `Genetic Variant`, fill = value)) +
  geom_tile()+
  scale_fill_manual(values=cols)+
  theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

func_binary <- function(x) {
  x = ifelse(x =='TRUE', 1, 0)
  return(x)
}

df_mutate[2:ncol(df_mutate)]<- lapply(df_mutate[2:ncol(df_mutate)], FUN = func_binary)

df_mutate2<-df_mutate[,-c(1)]
##########################################3
df_blast<-df%>%dplyr::select(c('Blast'))
df_blast[1]<- lapply(df_blast[1], FUN = function(y){as.numeric(y)})
df_blast[is.na(df_blast)] <- 0
df_blast<-data.matrix(df_blast)

df_AML<-df_Cell_Blast%>%dplyr::select(c(1))
df_AML[1]<- lapply(df_AML[1], FUN = function(y){as.numeric(y)})
df_AML[is.na(df_AML)] <- 0
df_AML<-data.matrix(df_AML)

df_LSC<-df_Cell_Blast%>%dplyr::select(c(2))
df_LSC[1]<- lapply(df_LSC[1], FUN = function(y){as.numeric(y)})
df_LSC[is.na(df_LSC)] <- 0
df_LSC<-data.matrix(df_LSC)

df_CC<-df_Cell_Blast%>%dplyr::select(c(3))
df_CC[1]<- lapply(df_CC[1], FUN = function(y){as.numeric(y)})
df_CC[is.na(df_CC)] <- 0
df_CC<-data.matrix(df_CC)


df_MRD<-df%>%dplyr::select(c("MRD"))
df_MRD[1]<- lapply(df_MRD[1], FUN = function(y){as.numeric(y)})
df_MRD[is.na(df_MRD)] <- 0
df_MRD<-data.matrix(df_MRD)

df_Age<-df%>%dplyr::select(c("Age"))
df_Age[1]<- lapply(df_Age[1], FUN = function(y){as.numeric(y)})
df_Age[is.na(df_Age)] <- 0
df_Age<-data.matrix(df_Age)
###############
df_clinic<-df%>%dplyr::select(c("FAB_2","Alive","Relapsed","Remission"))
df_clinic[is.na(df_clinic)] <- 'none'
colnames(df_clinic)
clinic_col <- list(Alive = c("1" = 'white', "0" = "#B2182B", "none" = "grey"),
            Relapsed= c("1" = '#B2182B', "0" ="#17becf", "none" = "grey"),
            Remission = c("1" = 'white', "0" = "#B2182B", "none" = "grey"),
            FAB=c('M1'='#d62728',
              'M2'='#aa40fc',
              "M3"= '#8c564b',
              "M4"='#e377c2',
              "M4Eo"='#b5bd61',
              "M5a"='#17becf',
              "M5b"='#aec7e8',
              "M6"='#ffbb78',
              "MPAL"='#98df8a',
              "ND"='grey',
              "none"='grey'))
C_ha = HeatmapAnnotation(which = 'row',
                         col = clinic_col,
                         FAB=df_clinic[[1]],
                         Alive = df_clinic[[2]],
                         Relapsed = df_clinic[[3]],
                         Remission = df_clinic[[4]],
                         border = TRUE
)

mat<-matrix(nrow = nrow(df), ncol = 0)
rownames(mat)<-df$Sample

library(dendextend)
row_dend = as.dendrogram(hclust(as.dist(1-rows.cor)))
row_dend = color_branches(row_dend, k = 4) # `color_branches()` returns a dendrogram object

Heatmap(mat, name = "mat", cluster_rows = row_dend)

ht1=Heatmap(mat,
        cluster_rows = row_dend,
        row_dend_width = unit(4, "cm"),
        right_annotation = C_ha)
ht1
col_fun = colorRamp2(c(0, 1), c("white",'#B2182B'))

ht2=Heatmap(df_mutate2,
            show_column_dend = FALSE,
            cluster_columns=FALSE,
            cluster_rows = hclust.row,
            col=col_fun,
            rect_gp = gpar(col = "grey", lwd = 1),
            border = TRUE,
            width = ncol(df_mutate2)*unit(5, "mm"))+
  rowAnnotation(`Blast(%)` = anno_barplot(df_blast,add_numbers = TRUE,width = unit(4, "cm")))+
  #rowAnnotation(`AML(%)` = anno_barplot(df_AML,add_numbers = TRUE,width = unit(4, "cm")))+
  #rowAnnotation(`AML_CD34(%)` = anno_barplot(df_LSC,add_numbers = TRUE,width = unit(4, "cm")))+
  #rowAnnotation(`AML_CellCycle(%)` = anno_barplot(df_CC,add_numbers = TRUE,width = unit(4, "cm")))+
  rowAnnotation(`MRD(%)` = anno_barplot(df_MRD,add_numbers = TRUE,width = unit(4, "cm")))+
  rowAnnotation(`Age(y)` = anno_points(df_Age),width = unit(4, "cm"))

ht=ht1+ht2

draw(ht, auto_adjust = FALSE)
decorate_annotation("Age(y)", {
  grid.lines(unit(c(1, 1), "native"), c(1, 0),gp = gpar(col = "red"))
})
decorate_annotation("Age(y)", {
  grid.lines(unit(c(10, 10), "native"), c(1, 0),gp = gpar(col = "red"))
})

