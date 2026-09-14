# ====================================================================
# Scanpy_in_R.R  |  CD96 figure pipeline component
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript Scanpy_in_R.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

rm(list=ls())

library(Seurat)
library(reticulate)
#remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
# install SeuratDisk from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/seurat-disk', ref = 'develop')
library(scater)
library(tidyverse)
library(cowplot)
library(harmony)
library(ggrepel)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(MetBrewer)
# using the cowplot theme for ggplot
library(patchwork)
library(proxy)
library(corrplot)

library(presto)
#devtools::install_github("immunogenomics/presto")
library(ComplexHeatmap)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ComplexHeatmap")
library(circlize)


# enable parallel processing for network analysis (optional)
enableWGCNAThreads(nThreads = 8)
allowWGCNAThreads(nThreads = 8)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

setwd('/oak/stanford/groups/cgawad/home/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/')

use_condaenv(condaenv = "pyscenic_3.8", conda = "/oak/stanford/groups/cgawad/Sequencing_Analysis_Tools/scRNA_py3.7/bin/conda")

pd <- import("pandas")
an <- import("anndata")

adata<-an$read_h5ad('/oak/stanford/groups/cgawad/home/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/ALSF_AML_Combo_3500_with_raw_count.h5ad')
metadata<-as.data.frame(adata$obs)
embedding <- adata$obsm["X_umap"]
rownames(embedding) <- adata$obs_names$to_list()
colnames(embedding) <- c("umap_1", "umap_2")

rm(adata)

Raw_data<-Read10X(data.dir='/oak/stanford/groups/cgawad/home/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/matrix_seurat/AML_counts')

#test<-as.data.frame(Raw_data)
#unique(test[,1138])
#rm(test)

metadata$cell_id<-rownames(metadata)

sobj<-CreateSeuratObject(counts = Raw_data,meta.data = metadata)

sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)

sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sobj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sobj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
plot2


all.genes <- rownames(sobj)
sobj <- ScaleData(sobj, features = all.genes)
#install.packages("irlba", type = "source")
#remotes::install_version("Matrix", version = "1.6-4")

sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))
sobj <- RunHarmony(sobj, group.by.vars = c("Sample"),  reduction = "pca",  reduction.save = "harmony")
sobj <- RunUMAP(sobj, reduction = "pca", dims = 1:30)
sobj <- FindNeighbors(object = sobj, reduction = "pca")

#sobj.list <- SplitObject(sobj, split.by = "SampleID")


#sobj.list <- lapply(X = sobj.list, FUN = function(x) {
#  x <- NormalizeData(x)
#  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#})


#features <- SelectIntegrationFeatures(object.list = sobj.list)

#sobj.list <- lapply(X = sobj.list, FUN = function(x) {
#  x <- ScaleData(x, features = features, verbose = FALSE)
 # x <- RunPCA(x, features = features, verbose = FALSE)
#})

#anchors <- FindIntegrationAnchors(object.list = sobj.list, reference = c(1,2), reduction = "rpca",  dims = 1:50)
#sobj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

#sobj<-sobj.integrated

#DefaultAssay(sobj) <- "integrated"

sobj<- ScaleData(sobj)
sobj<- RunPCA(sobj, features = features)
sobj<- RunUMAP(sobj, reduction = "pca", dims = 1:30)
sobj<- FindNeighbors(sobj, reduction = "pca", dims = 1:30)
sobj<- FindClusters(sobj, resolution = 0.5)

p1 <- DimPlot(sobj, reduction = "umap", group.by = "SampleID")
p2 <- DimPlot(sobj, reduction = "umap", group.by = "Cell_Type", label = TRUE,repel = TRUE)
p1
p2

#You can increase the strength of alignment by increasing the k.anchor parameter,
#which is set to 5 by default. 
#Increasing this parameter to 20 will assist in aligning these populations.

sobj <- SetIdent(sobj, value = "leiden")

sobj@reductions$umap<- CreateDimReducObject(embedding=embedding, key = "umap_",assay = 'RNA')


DimPlot(sobj, reduction = "umap")

p <- DimPlot(sobj, group.by='PAC_anno') +
  umap_theme() + ggtitle('PAC Control Cortex') 

p


saveRDS(sobj, "/oak/stanford/groups/cgawad/home/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/ALSF_AML_Leukemia_NormalHSPC_integrated.rds")
###Assign scanpy umap to seurat object
#sobj <- readRDS("/oak/stanford/groups/cgawad/home/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/ALSF_AML_Combo_integrated.rds")

  
sobj <- readRDS("/oak/stanford/groups/cgawad/home/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/ALSF_AML_Leukemia_NormalHSPC_integrated.rds")

seurat_obj <- SetupForWGCNA(
  sobj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.005, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "LSC" # the name of the hdWGCNA experiment
)

rm(sobj)

seurat_obj$all_cells <- 'all'
# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("all_cells","Sample"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'all_cells', # set the Idents of the metacell seurat object
)


#seurat_obj <- MetacellsByGroups(
 # seurat_obj = seurat_obj,
  #group.by = c('LSC_Cell_Type_2',"Sample"), # specify the columns in seurat_obj@meta.data to group by
  #k = 5, # nearest-neighbors parameter
  #min_cells = 10,
  #max_shared = 10, # maximum number of shared cells between two metacells
  #ident.group = 'LSC_Cell_Type_2' # set the Idents of the metacell seurat object
#)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "all", # the name of the group of interest in the group.by column
  group.by="all_cells",
  assay= 'RNA'# the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
)




# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj
  #setDatExpr = FALSE# set this to FALSE since we did this above
)


saveRDS(seurat_obj,"/oak/stanford/groups/cgawad/home/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/ALSF_AML_Leukemia_NormalHSPC_integrated_hdWGCNA.rds")

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(seurat_obj)


#Construct co-expression network
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=7,
  setDatExpr=FALSE,
  tom_name = 'AML' ,
  overwrite_tom = TRUE# name of the topoligical overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='AML hdWGCNA Dendrogram')
saveRDS(seurat_obj,"/oak/stanford/groups/cgawad/home/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/ALSF_AML_Leukemia_NormalHSPC_integrated_hdWGCNA.rds")
sobj <- readRDS("/oak/stanford/groups/cgawad/home/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/ALSF_AML_Leukemia_NormalHSPC_integrated_hdWGCNA.rds")
#sobj@meta.data<-metadata

##############################################
# get the module table
modules <- GetModules(seurat_obj)
mods <- unique(modules$module)

# make a table of the module-color pairings
mod_colors_df <- dplyr::select(modules, c(module, color)) %>%
  distinct %>% arrange(module)
rownames(mod_colors_df) <- mod_colors_df$module

# print the dataframe
mod_colors_df

# load MetBrewer color scheme pakckage
library(MetBrewer)

# get a table of just the module and it's unique color
mod_color_df <- GetModules(seurat_obj) %>%
  dplyr::select(c(module, color)) %>%
  distinct %>% arrange(module)

# the number of unique modules (subtract 1 because the grey module stays grey):
n_mods <- nrow(mod_color_df) - 1

# using the "Signac" palette from metbrewer, selecting for the number of modules
new_colors <- paste0(met.brewer("Signac", n=n_mods))

# reset the module colors
seurat_obj <- ResetModuleColors(seurat_obj, new_colors)
#################################################
PlotDendrogram(seurat_obj, main='AML hdWGCNA Dendrogram')
# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="SampleID"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)



# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group_name = "all",
  group.by = "all_cells"
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Module"
)


p <- PlotKMEs(seurat_obj, 
              n_hubs = 10,
              text_size = 3,
              ncol=6)

p


# get the module assignment table:
modules <- GetModules(seurat_obj)

# show the first 6 columns:
head(modules[,1:6])


# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 200)

head(hub_df)

write_csv2(hub_df, 'AML_Module_hub_gene_list.csv')

saveRDS(seurat_obj,
  "/oak/stanford/groups/cgawad/home/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/ALSF_AML_Leukemia_NormalHSPC_integrated_hdWGCNA.rds")
##########################################################################################
# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)

# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
#library(UCell)
# seurat_obj <- ModuleExprScore(
#   seurat_obj,
#   n_genes = 25,
#   method='UCell'
# )

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE,# order so the points with highest hMEs are on top
  #ucell = TRUE
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hub gene scores
  #order=TRUE,
  order='shuffle', # order so cells are shuffled
  #ucell = TRUE # depending on Seurat vs UCell for gene scoring
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

#############################
# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

PAC_object<-subset(x=seurat_obj,subset=PAC_anno%in%c("0_HSPC","Myeloid_Pro",'FPAC_1', 'FPAC_2', 'FPAC_3',
                                                     'FPAC_4',
                                                     'PPAC_1', 'PPAC_2', 'PPAC_3', 'PPAC_4', 'PPAC_5'
                                              ))



seurat_obj@meta.data$Prognosis

# plot with Seurat's DotPlot function

p <- DotPlot(seurat_obj,
             features=mods,
             group.by ='PAC_anno')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  xlab('')+ylab('')+
  theme(axis.text.x = element_text( size=9, hjust=0.8,angle = 90),
        axis.text.y = element_text(#face="bold",
                                   size=9
                                   ),
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9))+
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output
p

saveRDS(p, file = "PAC_anno_Module_data.rds")

p <- VlnPlot(
  seurat_obj,
  features = 'Module10',
  group.by = 'Cytogenetic',
  pt.size = 0 # don't show actual data points
)

# add box-and-whisker plots on top:
p <- p + geom_boxplot(width=.25, fill='white')

# change axis labels and remove legend:
p <- p + xlab('') + ylab('hME') + NoLegend()

# plot output
p
##############################################
library('igraph')
ModuleNetworkPlot(seurat_obj,
                  vertex.label.cex = 0.6,
                  vertex.size = 6)
# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all',
  vertex.label.cex = 0.3,
  hub.vertex.size = 4
)

dev.off()

#g <- HubGeneNetworkPlot(seurat_obj, return_graph=TRUE)

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 50, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.2 # min distance between points in UMAP space
)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

dev.off() 

options(future.globals.maxSize = 8000* 1024^3)


ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=FALSE,# TRUE if we sample edges randomly, FALSE if we take the top edges
  edge_prop=0.2, # proportion of edges to sample (20% here)
  label_hubs=5 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE,
  vertex.label.cex = 0.08,
  )


#g <- ModuleUMAPPlot(seurat_obj,  return_graph=TRUE)
df<-seurat_obj@meta.data
colnames(df)

df<-df%>%subset(select=c("cell_id","PAC_anno","Cytogenetic",
                         "Module10" ,                    
                         "Module4",                       "Module7" ,                      "Module17" ,                    
                         "Module8",                      "Module12",                      "Module18",                     
                         "Module9",                       "Module15",                      "Module5",                      
                         "Module3",                       "Module2",                       "Module1",                      
                         "Module14",                      "Module16",                     
                         "Module11",                      "Module6",                       "Module13"  ))
write_tsv(df, "ALSF_AML_module_Cytogenetic_MEs.tsv")
###################################################

seurat_obj$Prognosis<- as.factor(seurat_obj$Prognosis)
seurat_obj$LSC_EPPERT_Score <- as.numeric(seurat_obj$LSC_EPPERT_Score)
seurat_obj$HSC_Score <- as.numeric(seurat_obj$HSC_Score)
seurat_obj$LSC17_Score<- as.numeric(seurat_obj$LSC17_Score)
seurat_obj$Prognosis_Score<- as.numeric(seurat_obj$Prognosis_Score)
seurat_obj$Relapse_Score <- as.numeric(seurat_obj$Relapse_Score)

cur_traits <- c("LSC_EPPERT_Score","HSC_Score","LSC17_Score")

cur_traits <- c("Prognosis_Score","Relapse_Score")


seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits,
  group.by='PAC_anno'
)


PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = 'yellow',
  mid_color = 'black',
  low_color = 'purple',
  plot_max = 0.2,
  combine=TRUE
)
###################################################


#g <- ModuleUMAPPlot(seurat_obj,  return_graph=TRUE)
df<-MEs
df$Sample<-seurat_obj@meta.data$PAC_anno
#df$Cell_Type<-seurat_obj@meta.data$Cell_Type
#df<-filter(df, Cell_Type=='AML')
colnames(df)
write_tsv(df, "//oak/stanford/groups/cgawad/home/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/H5AD/ALSF_AML_18_module_MEs_by_PAC_anno.tsv")
#####expression mean######
df_mean<-df%>%group_by(Sample)%>%
  summarise(across(`Module1`:`Module18`,~ mean(.x, na.rm = TRUE)))

write_tsv(df_mean, "ALSF_AML_module_mean.tsv")
#####expression fraction######
df_counts<-df%>%group_by(Sample)%>%
  summarise(across(`Module1`:`Module18`,  ~ sum(.x > 0, na.rm = TRUE)))

write_tsv(df_counts, "ALSF_AML_module_sum.tsv")


df_counts<-df_counts%>%mutate(total = rowSums(across(where(is.numeric))))

df_pct<-df_counts%>%summarise(across(`Module1`:`Module18`, ~.x/total))
df_pct$Sample<-df_counts$Sample
write_tsv(df_pct, "ALSF_AML_module_pct.tsv")
##############################################
AML_object<-subset(x=seurat_obj,subset=Cell_Type%in%c("AML"))
p <- DotPlot(AML_object, features=mods, group.by = 'Sample')
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')
p
######################################################
# gene enrichment packages
library(enrichR)
library(GeneOverlap)

dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# perform enrichment tests
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)
# make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)

EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=3,# number of terms for each module
)
########################################################
mods<-c('AML-Module5','AML-Module7','AML-Module10','AML-Module17')
mods<-c('AML-Module6','AML-Module14','AML-Module13','AML-Module9')

plot_df <- enrich_df %>%
  subset(db =="GO_Biological_Process_2021" & module %in% mods) %>%
  group_by(module) %>%
  top_n(10, wt=Combined.Score)

plot_df<-merge(plot_df, mod_color_df)

plot_df$module <- factor(
  as.character(plot_df$module),
  levels=levels(modules$module)
)
plot_df <- arrange(plot_df, module)
plot_df$Term <- factor(
  as.character(plot_df$Term),
  levels = unique(as.character(plot_df$Term))
)

plot_df$Combined.Score <- log(plot_df$Combined.Score)

lab <- 'Enrichment\nlog(combined score)'


p <- plot_df  %>%
  ggplot(aes(x=module, y=Term)) +
  geom_point(aes(size=Combined.Score), color=plot_df$color) +
  RotatedAxis() +
  ylab('') + xlab('') +
  theme_minimal_grid()+
  labs(size=lab) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))+
  #ggtitle("AML Module 6 9 13 14") +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.text.x = element_text( size=9, hjust=0.8,angle = 90))



p


#########################################################
ModuleCorrelogram(seurat_obj)
#DME analysis comparing two groups
HBM<-c("0_HealthyBM1", "0_HealthyBM2")

unique(seurat_obj@meta.data[,c('Prognosis','SampleID')])

unique(seurat_obj@meta.data$PAC_anno)

`%ni%` <- Negate(`%in%`)

group1 <- seurat_obj@meta.data %>% subset(Prognosis == 'Alive') %>% rownames
group2 <- seurat_obj@meta.data %>% subset(Prognosis == 'Deceased')%>% rownames

AML_object<-subset(x=seurat_obj,subset=lineage%in%c('AML'))

group1 <- AML_object@meta.data %>% subset(Prognosis == 'Alive') %>% rownames
group2 <- AML_object@meta.data %>% subset(Prognosis == 'Deceased')%>% rownames

head(group1)

DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name='LSC'
)

head(DMEs)
DMEs$anno <- ifelse(DMEs$p_val_adj < 0.05, DMEs$module, '')

modules <- GetModules(seurat_obj, "LSC") %>% subset(module != 'grey') %>% mutate(module=droplevels(module))
module_colors <- modules %>% dplyr::select(c(module, color)) %>% distinct

# module names
mods <- levels(modules$module)
mods <- mods[mods %in% DMEs$module]
mod_colors <- module_colors$color; names(mod_colors) <- as.character(module_colors$module)

p1 <- ggplot(DMEs, aes(avg_log2FC, -log(p_val_adj,10), fill=module, color=module)) + # -log10 conversion  
  geom_point(size=5, pch=21, color='black',) +
  scale_fill_manual(values=mod_colors) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))+
  theme(
    panel.border = element_rect(color='black', fill=NA, size=1),
    panel.grid.major = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position='bottom'
  ) + NoLegend()

p1 <- p1 +
  geom_vline(xintercept=0, linetype='dashed', color='grey75', alpha=0.8) +
  geom_rect(
    data=DMEs[1,],
    aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=-log10(0.05)), fill='grey75', alpha=0.8, color=NA)+
  geom_text_repel(aes(label=anno), color='black', min.segment.length=0, max.overlaps=Inf, size=3)


p1





PlotDMEsVolcano(
  DMEs,
  mod_point_size = 3,
  label_size = 3,
  wgcna_name='LSC'
)

AML_object<-subset(x=seurat_obj,
                   subset=PAC_anno%in%c('0_HSPC','Myeloid_Pro','AML','AML-MKI67','AML-PCNA',
                                        "FPAC_1","FPAC_2",
                                        "PPAC_1","PPAC_2","PPAC_3","PPAC_4","PPAC_5"))

control_object<-subset(x=seurat_obj,
                       subset=Sample%in%HBM & subset=PAC_anno%in%c('0_HSPC','Myeloid_Pro','CD14_Monocyte'))

#AML_object<-subset(x=seurat_obj,subset=Sample%ni%HBM)

DMEs_all <- FindAllDMEs(
  AML_object,
  group.by = 'PAC_anno',
  wgcna_name = 'LSC'
)

head(DMEs_all)

DMEs_all<-DMEs_all%>%mutate(p_val_adj=DMEs_all$p_val_adj+1e-200)

head(DMEs_all)
  
p <- PlotDMEsVolcano(
  DMEs_all,
  wgcna_name = 'LSC',
  mod_point_size = 2,
  label_size = 2,
  # plot_labels=FALSE,
  # show_cutoff=FALSE
)

DMEs_all$anno <- ifelse(DMEs_all$p_val_adj < 0.05, DMEs_all$module, '')

modules <- GetModules(seurat_obj, "LSC") %>% subset(module != 'grey') %>% mutate(module=droplevels(module))
module_colors <- modules %>% dplyr::select(c(module, color)) %>% distinct

# module names
mods <- levels(modules$module)
mods <- mods[mods %in% DMEs_all$module]
mod_colors <- module_colors$color; names(mod_colors) <- as.character(module_colors$module)



p1 <- ggplot(DMEs_all, aes(avg_log2FC, -log(p_val_adj,10), fill=module, color=module)) + # -log10 conversion  
  geom_point(size=4, pch=21, color='black',) +
  scale_fill_manual(values=mod_colors) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))+
  theme(
    panel.border = element_rect(color='black', fill=NA, size=1),
    panel.grid.major = element_blank(),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.position='bottom'
  ) + NoLegend()

p1 <- p1 +
  geom_vline(xintercept=0, linetype='dashed', color='grey75', alpha=0.8) +
  geom_rect(
    data=DMEs[1,],
    aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=-log10(0.05)), fill='grey75', alpha=0.8, color=NA)+
  geom_text_repel(aes(label=anno), color='black', min.segment.length=0, max.overlaps=Inf, size=3)


p1+facet_wrap(~forcats::fct_relevel(group, '0_HSPC','Myeloid_Pro','AML','AML-MKI67','AML-PCNA',
                                  "FPAC_1","FPAC_2",
                                  "PPAC_1","PPAC_2","PPAC_3","PPAC_4","PPAC_5"),
                      ncol=4,scales='free_y')
# facet wrap by each cell type
p + facet_wrap(~group, ncol=3)

###########################################################
library(JASPAR2020)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)


# get the pfm from JASPAR2020 using TFBSTools
pfm_core <- TFBSTools::getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# run the motif scan with these settings for the mouse dataset
seurat_obj <- MotifScan(
  seurat_obj,
  species_genome = 'hg38',
  pfm = pfm_core,
  EnsDb = EnsDb.Hsapiens.v86
)
dim(GetMotifMatrix(seurat_obj))


# TF target genes
target_genes <- GetMotifTargets(seurat_obj)

# check target genes for one TF:
head(target_genes$SOX9)


seurat_obj<- OverlapModulesMotifs(seurat_obj)

# look at the overlap data
head(GetMotifOverlap(seurat_obj))

# plot the top TFs overlapping with
MotifOverlapBarPlot(
  seurat_obj,
  #motif_font = 'xkcd_regular',
  outdir = 'motifs/MotifOverlaps/',
  plot_size=c(5,6)
)


library(UCell)
seurat_obj <- MotifTargetScore(
  seurat_obj,
  method='UCell'
)

df <- GetMotifOverlap(seurat_obj)

cur_df <- df %>% subset(tf == 'SOX9')

plot_var <- 'odds_ratio'
p <- cur_df %>%
  ggplot(aes(y=reorder(module, odds_ratio), x=odds_ratio)) +
  geom_bar(stat='identity', fill=cur_df$color) +
  geom_vline(xintercept = 1, linetype='dashed', color='gray') +
  geom_text(aes(label=Significance), color='black', size=3.5, hjust='center') +
  ylab('') +
  xlab("Odds Ratio") +
  ggtitle("SOX9 overlap") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


png(paste0(fig_dir, 'Sox9_motif_overlap_or.png'), width=3, height=4, units='in', res=400)
p
dev.off()


ModuleTFNetwork(
  seurat_obj,
  edge.alpha=0.75,
  cor_thresh = 0.75,
  tf_name = "OLIG1",
  tf_gene_name = "OLIG1",
  tf_x = -7,
  tf_y = -7
)


