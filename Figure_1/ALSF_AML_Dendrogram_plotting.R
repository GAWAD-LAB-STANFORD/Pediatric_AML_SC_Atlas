# ====================================================================
# ALSF_AML_Dendrogram_plotting.R  |  Figure 1E (per-sample dendrogram)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : AML-sample-infor.csv, AML_scRNA_Cell_summary_Cell_Type.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript ALSF_AML_Dendrogram_plotting.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================


rm(list=ls())
suppressPackageStartupMessages({
  library(tidyverse)
})
library(ggplot2)
library(reshape2)
library(dplyr)
detach(package:plyr)

setwd("~/Downloads/ALSF_AML_plot")
df_pca <- read.table("pca_Combo_AML_with_normal_HSPC_by_lineage_Sample_for_dendrogram_rePCA.csv",sep = ",", header = TRUE)
df_pca[2:ncol(df_pca)]<- lapply(df_pca[2:ncol(df_pca)], FUN = function(y){as.numeric(y)})
df_pca<-df_pca%>% separate(col=1, sep = "-", c("Sample", "Drop1"))
df_pca[3:ncol(df_pca)]<- lapply(df_pca[3:ncol(df_pca)], FUN = function(y){as.numeric(y)})
################################################################################################
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


df <-merge(df_clic,df_w, by="Sample", all = T)
df <-merge(df_clic,df_Cell_mean[,-2], by="Sample", all = T)
df <-merge(df_clic,df_RSS[,-2], by="Sample", all = T)
rownames(df)<-df$Sample

df_clic_subset<-df_clic[,c('Sample','Cytogenetic')]
###################################################################
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


rownames(df_w) <-df_w$Sample
data<-df_w[,-c(1)]
data[is.na(data)]<-"none"
data.dist <- vegdist(data[,c(2:ncol(data))],method = "jaccard")
row.clus <- hclust(data.dist,'ward.D2')

dhc <- as.dendrogram(row.clus)
plot(row.clus,hang = -1, cex = 0.6)
#change the attribute of each leaf. This can be done using the dendrapply function. 
#create a function that # # add 3 attributes to the leaf : one for the color (“lab.col”) ,
#one for the font “lab.font” and one for the size (“lab.cex”)

rownames(df_pca) <-df_pca$Sample
Sample_list<-c( "AML335","AML647","AML882","AML948","AML1355",      
              "AML3082","AML3121","AML3210" ,"AML3371" , "AML3492" , "AML4000",      
              "AML4010","AML4035","AML4062" , "AML4068" ,"AML4090",      
              "AML4102", "AML4116","AML4127","AML4192" , "AML4226", "AML4232",      
              "AML4239","AML4264","AML4271" , "AML4304", "AML4363" )
df_pca<-filter(df_pca, Sample%in%Sample_list)

df <-merge(df_clic_subset,df_w, by="Sample")
rownames(df)<-df$Sample
data<-df[,-c(1,2)]
rownames(data) <-df$Sample
data[is.na(data)]<-"none"


data.dist <- vegdist(data[,c(2:ncol(data))],method = "jaccard")
row.clus <- hclust(data.dist,'ward.D2')

unique(df_clic$Cytogenetic)

dhc <- as.dendrogram(row.clus)
plot(row.clus,hang = -1, cex = 0.6)

ordered_labels <- labels(dhc)[order.dendrogram(dhc)]

print(ordered_labels)

color<-c('#d62728',
'#aa40fc',
'#8c564b',
'#e377c2',
'#b5bd61',
'#17becf',
'#aec7e8',
'#ffbb78',
'#98df8a',
'#ff9896',
'darkblue',
"grey",
"black")

col_FAB=""
col_Prognosis=""
col_Relapsed=""
col_Cytogenetic=""

colLab<<-function(n){
  if(is.leaf(n)){
    
    #I take the current attributes
    a=attributes(n)
    
    #I deduce the line in the original data, and so the treatment and the specie.
    ligne=match(attributes(n)$label,df[,1])

    Cytogenetic=df[ligne,2];
      if(Cytogenetic=="PML/RARA"){col_Cytogenetic='#D62728'}
      if(Cytogenetic=="CBFB/MYH11"){col_Cytogenetic='#AA40FC'}
      if(Cytogenetic=="t(2;3)(p15;q26.2)"){col_Cytogenetic='#8C564B'}
      if(Cytogenetic=="MLLr"){col_Cytogenetic='#E377C2'}
      if(Cytogenetic=="CN"){col_Cytogenetic='#B5BD61'}
      if(Cytogenetic=="Tri(15)"){col_Cytogenetic='#17BECF'}
      if(Cytogenetic=="MYB/GATA1"){col_Cytogenetic='#AEC7E8'}
      if(Cytogenetic=="Tri(8)/MLLr"){col_Cytogenetic='#FFBB78'}
      if(Cytogenetic=="RUNX1/RUNX1T1"){col_Cytogenetic='#98DF8A'}
      if(Cytogenetic=="NUP98/NSD1"){col_Cytogenetic="black"}
      if(Cytogenetic=="Tri(8)"){col_Cytogenetic="grey"}
      if(Cytogenetic=='del7q'){col_Cytogenetic='#279E68'}
      if(Cytogenetic=='BCR/ABL'){col_Cytogenetic='#FF9896'}
      if(Cytogenetic=='t(7;14)(q21;q32)'){col_Cytogenetic='darkblue'}
    
    
   # Prognosis=df[ligne,30];
    #  if(Prognosis=='0'){col_Prognosis='#279E68'}
    #  if(Prognosis=='1'){col_Prognosis='#FF9896'}
    #  if(Prognosis=='none'){col_Prognosis="black"}
    
   # Relapsed=df[ligne,25];
   # if(Relapsed=='0'){col_Relapsed='#279E68'}
   # if(Relapsed=='1'){col_Relapsed='#FF9896'}
   # if(Relapsed=='none'){col_Relapsed="black"}
    
    #Modification of leaf attribute
    attr(n,"nodePar")<-c(a$nodePar,list(cex=1.5,lab.cex=0.8,pch=20,
                                        col=col_Cytogenetic,
                                        lab.col=col_Cytogenetic,
                                        lab.font=1,lab.cex=1))
  }
  return(n)
}


dL <- dendrapply(dhc, colLab)

plot(dL,main="structure of the population")


legend(28, 45,
       legend = c("del7q",'BCR/ABL',"PML/RARA","t(2;3)(p15;q26.2)","CBFB/MYH11",
                  "MLLr","CN","Tri(15)", "MYB/GATA1", "Tri(8)/MLLr",  "RUNX1/RUNX1T1" ,
                  "NUP98/NSD1","Tri(8)","t(7;14)(q21;q32)"), 
       col = c('#279E68', '#FF9896' , '#D62728' ,'#AA40FC' ,'#8C564B','#E377C2',
               '#B5BD61','#17BECF','#AEC7E8','#FFBB78','#98DF8A','black','grey',
               'darkblue'), 
       pch = c(8,8,rep(20,10)),
       bty = "n",  
       pt.cex = 1.5, 
       cex = 0.6 , 
       text.col = "black", horiz = FALSE, 
       inset = c(0, 0.1) # Distance from the margin as a fraction of the plot region
       )
####################################################################
mydata<-df_pca[3:52]
rownames(mydata)<-df_pca$Sample
rows.cor <- cor(t(mydata), use = "complete.obs", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor))
hclust.row$labels
########################################
mydata<-df_Cell_mean[3:235]
rownames(mydata)<-df_Cell_mean$Sample
rows.cor <- cor(t(mydata), use = "complete.obs", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor))
###############################################
mydata<-df_RSS[-c(1,2)]
rownames(mydata)<-df_RSS$Sample
rows.cor <- cor(t(mydata), use = "complete.obs", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor))
##################################################################
mydata<-df_w[2:ncol(df_w)]
rownames(mydata)<-df_w$Sample
rows.cor <- cor(t(mydata), use = "complete.obs", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor))
####################
colnames(df)
df_cytoge<-df%>%dplyr::select(c("Sample", "Cytogenetic"))
df_cytoge<-df_cytoge%>%mutate(value='TRUE')
df_cytoge_w<-spread(df_cytoge, key='Cytogenetic', value='value',fill='none')
df[6:30]
df_mutate<-cbind(df_cytoge_w[1:15],df[6:30])
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

df_AML<-df%>%dplyr::select(c(1))
df_AML[1]<- lapply(df_AML[1], FUN = function(y){as.numeric(y)})
df_AML[is.na(df_AML)] <- 0
df_AML<-data.matrix(df_AML)

df_LSC<-df%>%dplyr::select(c(2))
df_LSC[1]<- lapply(df_LSC[1], FUN = function(y){as.numeric(y)})
df_LSC[is.na(df_LSC)] <- 0
df_LSC<-data.matrix(df_LSC)

df_CC<-df%>%dplyr::select(c(3))
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
colnames(df_clinic)
clinic_col <- list(Alive = c("1" = 'white', "0" = "#B2182B", "none" = "grey"),
            Relapsed= c("1" = 'white', "0" ="#98df8a", "none" = "grey"),
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

ht1=Heatmap(mat,
        cluster_rows = hclust.row,
        row_dend_width = unit(4, "cm"),
        right_annotation = C_ha)
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
################################

ht=ht1+ht2




ht = Heatmap(heat.in,cluster_columns = col_dend,
             cluster_rows = row_dend,width = ncol(heat.in)*unit(3, "mm"), 
             height = nrow(heat.in)*unit(1.5, "mm"),col=col_fun,
             column_names_gp = grid::gpar(fontsize = 6),
             row_names_gp = grid::gpar(fontsize = 4), 
             rect_gp = gpar(col = "#EEEEEE", lwd = 0.5),
             heatmap_legend_param = list(color_bar = "discrete",
                                         title = "Genotype", 
                                         at = c(0, 1, 2, 4), 
                                         labels = gt_render(c("Uncovered","WT","Het","Homo")),
                                         legend_gp = gpar(fill = col_fun)
             ),
             bottom_annotation = Bulk_ha,
             left_annotation = CNV_ha)
draw(ht, show_annotation_legend = FALSE, annotation_legend_side = "bottom")
dev.off()

del7_gene<-
