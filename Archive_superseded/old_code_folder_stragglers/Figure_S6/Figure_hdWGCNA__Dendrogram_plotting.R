# ====================================================================
# Figure_hdWGCNA__Dendrogram_plotting.R  |  Figures S6/S7 (module dendrogram / dot plot)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript Figure_hdWGCNA__Dendrogram_plotting.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
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
library(MetBrewer)

#For tree Drawing
require(phangorn)
require(ape)
require(adegenet)
library(dplyr)
library(tibble)
library(dendextend)
library(ComplexHeatmap)
#BiocManager::install("ComplexHeatmap")
library(circlize)

detach(package:plyr)

setwd("~/Downloads/")

PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

df_pca <- read.table("/Users/yakun/Downloads/ALSF_AML_module_Cytogenetic_MEs.tsv",sep = "\t", header = TRUE)
df_pca[4:ncol(df_pca)]<- lapply(df_pca[4:ncol(df_pca)], FUN = function(y){as.numeric(y)})

df_ME<-df_pca%>%group_by(PAC_anno)%>%
  mutate(across(starts_with("Module"), mean, .names="mean_{col}"))

df_Pct<-df_pca%>%group_by(PAC_anno)%>%
  mutate(across(starts_with("Module"), ~PercentAbove(.x,0), .names="pct_{col}"))

df_ME_rd<-unique(df_ME[,c(2,22:39)])
df_Pct_rd<-unique(df_Pct[,c(2,22:39)])

df_ME_T<-as.data.frame(t(df_ME_rd))
colnames(df_ME_T)<-df_ME_T[1,]
df_ME_T<-df_ME_T[-1,]
#######################################

mydata<-df_ME_T
mydata[1:ncol(mydata)]<- lapply(mydata[1:ncol(mydata)], 
                                FUN = function(y){as.numeric(y)})

rows.cor <- cor(t(mydata), use = "everything", method = "pearson")
hclust.row <- hclust(as.dist(1-rows.cor))
hclust.row$labels

mat<-matrix(nrow = nrow(mydata), ncol = 0)
rownames(mat)<-rownames(mydata)

library(dendextend)
row_dend = as.dendrogram(hclust(as.dist(1-rows.cor)))
row_dend = color_branches(row_dend, k = 5) # `color_branches()` returns a dendrogram object

Heatmap(mat, name = "mat", cluster_rows = row_dend)

ht1=Heatmap(mat,
            cluster_rows = hclust.row,
            row_dend_width = unit(4, "cm"))
ht1

colnames(df_pca)


df_ME_rd<-unique(df_ME[,c(2,22:39)])
df_Pct_rd<-unique(df_Pct[,c(2,22:39)])


library(Seurat)
library(tidyverse)
#install.packages("devtools")
#install.packages("remotes")
remotes::install_github("https://github.com/immunogenomics/presto.git")
library(presto)



P<-readRDS(file = "/Users/yakun/Downloads/PAC_anno_Module_data.rds")
P$data 

df<- P$data
head(df)

exp_mat<-df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()

head(exp_mat)

## the matrix for the percentage of cells express a gene

percent_mat<-df %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

head(percent_mat)


library(viridis)
library(Polychrome)

Polychrome::swatch(viridis(20))

## get an idea of the ranges of the matrix
quantile(exp_mat, c(0.1, 0.5, 0.9, 0.99))

## any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-0.5,0.5, 2), viridis(20)[c(1,10, 20)])


cell_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}

## also do a kmeans clustering for the genes with k = 4
Heatmap(exp_mat,
        heatmap_legend_param=list(title="expression"),
        column_title = "clustered dotplot", 
        col=col_fun,
        rect_gp = gpar(type = "none"),
        cell_fun = cell_fun,
        row_names_gp = gpar(fontsize = 5),
        row_km = 5,
        border = "black")
###############################
layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= pindex(percent_mat, i, j)/100 * unit(2, "mm"),
              gp = gpar(fill = col_fun(pindex(exp_mat, i, j)), col = NA))}


lgd_list = list(
  Legend( labels = c(0,0.25,0.5,0.75,1), title = "Percentage",
          nrow = 1,
          graphics = list(
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.25 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.5 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 0.75 * unit(2, "mm"),
                                             gp = gpar(fill = "black")),
            function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                             gp = gpar(fill = "black"))),
            
  ))

set.seed(123)    
hp<- Heatmap(exp_mat,
             heatmap_legend_param=list(title="Expression",
                                       legend_direction = "horizontal"),
             column_title = "", 
             clustering_distance_rows = "pearson",
             clustering_distance_columns = "spearman",
    # ("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall")
             clustering_method_rows = "complete",
             clustering_method_columns = "complete",
             col=col_fun,
             rect_gp = gpar(type = "none"),
             layer_fun = layer_fun,
             row_names_gp = gpar(fontsize = 10),
             row_km = 4,
             border = "black"   )

draw( hp, 
      heatmap_legend_side = "bottom", 
      annotation_legend_list = lgd_list, 
      annotation_legend_side="bottom")
