# ====================================================================
# RNA-seq.R  |  CD96 figure pipeline component
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : (reads upstream object / raw matrix — see body)
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript RNA-seq.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

library(dplyr)
library(Seurat)
library(patchwork)
library(devtools)
library(URD)


# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/oak/stanford/groups/cgawad/R_Testing/Seurat/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
## An object of class Seurat 
## 13714 features across 2700 samples within 1 assay 
## Active assay: RNA (13714 features, 0 variable features)
##What does data in a count matrix look like?

##Standard pre-processing workflow
##The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. 
##These represent the selection and filtration of cells based on QC metrics, data normalization and scaling, 
##and the detection of highly variable features.

##QC and selecting cells for further analysis
##Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. A 
##few QC metrics commonly used by the community include

##The number of unique genes detected in each cell.
##Low-quality cells or empty droplets will often have very few genes
##Cell doublets or multiplets may exhibit an aberrantly high gene count
##Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
##The percentage of reads that map to the mitochondrial genome
##Low-quality / dying cells often exhibit extensive mitochondrial contamination
##We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the 
##percentage of counts originating from a set of features
##We use the set of all genes starting with MT- as a set of mitochondrial genes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

##Where are QC metrics stored in Seurat?
##In the example below, we visualize QC metrics, and use these to filter cells.

##We filter cells that have unique feature counts over 2,500 or less than 200
##We filter cells that have >5% mitochondrial counts
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

##Normalizing the data
##After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in pbmc[["RNA"]]@data.

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
##For clarity, in this previous line of code (and in future commands), we provide the default values for certain parameters in the function call. However, this isn’t required and the same behavior can be achieved with:

pbmc <- NormalizeData(pbmc)
##Identification of highly variable features (feature selection)
##We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). We and others have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.

##Our procedure in Seurat is described in detail here, and improves on previous versions by directly modeling the mean-variance relationship inherent in single-cell data, and is implemented in the FindVariableFeatures() function. By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)



##Scaling the data
##Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:

##Shifts the expression of each gene, so that the mean expression across cells is 0
##Scales the expression of each gene, so that the variance across cells is 1
##This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
##The results of this are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

##This step takes too long! Can I make it faster?
##  How can I remove unwanted sources of variation, as in Seurat v2?


##Perform linear dimensional reduction Next we perform PCA on the scaled data. By default, 
##only the previously determined variable features are used as input, but can be defined using 
##features argument if you wish to choose a different subset.

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

##Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
## PC_ 1 
## Positive:  CST3, TYROBP, LST1, AIF1, FTL 
## Negative:  MALAT1, LTB, IL32, IL7R, CD2 
## PC_ 2 
## Positive:  CD79A, MS4A1, TCL1A, HLA-DQA1, HLA-DQB1 
## Negative:  NKG7, PRF1, CST7, GZMB, GZMA 
## PC_ 3 
## Positive:  HLA-DQA1, CD79A, CD79B, HLA-DQB1, HLA-DPB1 
## Negative:  PPBP, PF4, SDPR, SPARC, GNG11 
## PC_ 4 
## Positive:  HLA-DQA1, CD79B, CD79A, MS4A1, HLA-DQB1 
## Negative:  VIM, IL7R, S100A6, IL32, S100A8 
## PC_ 5 
## Positive:  GZMB, NKG7, S100A8, FGFBP2, GNLY 
## Negative:  LTB, IL7R, CKB, VIM, MS4A7
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")


DimPlot(pbmc, reduction = "pca")


##In particular DimHeatmap() allows for easy exploration of the primary sources of heterogeneity 
##in a dataset, and can be useful when trying to decide which PCs to include for further downstream 
##analyses. Both cells and features are ordered according to their PCA scores. Setting cells to a number 
##plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. 
##Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated feature sets.

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)


DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

##Determine the ‘dimensionality’ of the dataset
##To overcome the extensive technical noise in any single feature for scRNA-seq data, 
##Seurat clusters cells based on their PCA scores, with each PC essentially representing a 
##‘metafeature’ that combines information across a correlated feature set. The top principal 
##components therefore represent a robust compression of the dataset. However, how many components 
##should we choose to include? 10? 20? 100?

##In Macosko et al, we implemented a resampling test inspired by the JackStraw procedure. We randomly 
##permute a subset of the data (1% by default) and rerun PCA, constructing a ‘null distribution’ of feature 
##scores, and repeat this procedure. We identify ‘significant’ PCs as those who have a strong enrichment of 
##low p-value features.

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)


##The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for 
##each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of features 
##with low p-values (solid curve above the dashed line). In this case it appears that there is a sharp drop-off in 
##significance after the first 10-12 PCs.

JackStrawPlot(pbmc, dims = 1:15)

##An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage 
##of variance explained by each one (ElbowPlot() function). In this example, we can observe an ‘elbow’ around PC9-10, 
##suggesting that the majority of true signal is captured in the first 10 PCs.

ElbowPlot(pbmc)


##Identifying the true dimensionality of a dataset – can be challenging/uncertain for the user. We therefore suggest these 
##three approaches to consider. The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, 
##and could be used in conjunction with GSEA for example. The second implements a statistical test based on a random null 
##model, but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that is 
##commonly used, and can be calculated instantly. In this example, all three approaches yielded similar results, but we might 
##have been justified in choosing anything between PC 7-12 as a cutoff.

##We chose 10 here, but encourage users to consider the following:

##Dendritic cell and NK aficionados may recognize that genes strongly associated with PCs 12 and 13 define rare immune subsets 
##(i.e. MZB1 is a marker for plasmacytoid DCs). However, these groups are so rare, they are difficult to distinguish from background 
##noise for a dataset of this size without prior knowledge. We encourage users to repeat downstream analyses with a different number 
##of PCs (10, 15, or even 50!). As you will observe, the results 
##often do not differ dramatically.
##We advise users to err on the higher side when choosing this parameter. For example, performing downstream analyses with only 5 PCs 
##does significantly and adversely affect results.

##Cluster the cells
##Seurat v3 applies a graph-based clustering approach, building upon initial strategies in (Macosko et al). Importantly, the distance metric 
##which drives the clustering analysis (based on previously identified PCs) remains the same. However, our approach to partitioning the 
##cellular distance matrix into clusters has dramatically improved. Our approach was heavily inspired by recent manuscripts which applied 
##graph-based clustering approaches to scRNA-seq data [SNN-Cliq, Xu and Su, Bioinformatics, 2015] and CyTOF data [PhenoGraph, Levine et al., 
##Cell, 2015]. Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between 
##cells with similar feature expression patterns, and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.

##As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between two 
##cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the FindNeighbors() 
##function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs).

##To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) or SLM [SLM, Blondel et al., 
##Journal of Statistical Mechanics], to iteratively group cells together, with the goal of optimizing the standard modularity function. 
##The FindClusters() function implements this procedure, and contains a resolution parameter that sets the ‘granularity’ of the downstream 
##clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.4-1.2 typically 
##returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters can 
##be found using the Idents() function.

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 2638
## Number of edges: 95927
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8728
## Number of communities: 9
## Elapsed time: 0 seconds
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
## AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1 AAACCGTGCTTCCG-1 
##                2                3                2                1 
## AAACCGTGTATGCG-1 
##                6 
## Levels: 0 1 2 3 4 5 6 7 8
##Run non-linear dimensional reduction (UMAP/tSNE)
##Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to 
##visualize and explore these datasets. The goal of these algorithms is to learn the underlying 
##manifold of the data in order to place similar cells together in low-dimensional space. Cells 
##within the graph-based clusters determined above should co-localize on these dimension reduction 
##plots. As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

##You can save the object at this point so that it can easily be loaded back in without 
##having to rerun the computationally intensive steps performed above, or easily shared 
##with collaborators.

saveRDS(pbmc, file = "/oak/stanford/groups/cgawad/R_Testing/Seurat/filtered_gene_bc_matrices/hg19/pbmc_tutorial.rds")
##Finding differentially expressed features (cluster biomarkers)
##Seurat can help you find markers that define clusters via differential expression. By default, it identifies positive and 
##negative markers of a single cluster (specified in ident.1), compared to all other cells. FindAllMarkers() automates this 
##process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

##The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, and 
##the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. 
##You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are 
##unlikely to be highly discriminatory. As another option to speed up these computations, max.cells.per.ident can be set. This will 
##downsample each identity class to have no more cells than whatever this is set to. While there is generally going to be a loss in 
##power, the speed increases can be significant and the most highly differentially expressed features will likely still rise to the top.

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
##             p_val avg_log2FC pct.1 pct.2    p_val_adj
## IL32 2.892340e-90  1.2013522 0.947 0.465 3.966555e-86
## LTB  1.060121e-86  1.2695776 0.981 0.643 1.453850e-82
## CD3D 8.794641e-71  0.9389621 0.922 0.432 1.206097e-66
## IL7R 3.516098e-68  1.1873213 0.750 0.326 4.821977e-64
## LDHB 1.642480e-67  0.8969774 0.954 0.614 2.252497e-63
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
##                       p_val avg_log2FC pct.1 pct.2     p_val_adj
## FCGR3A        8.246578e-205   4.261495 0.975 0.040 1.130936e-200
## IFITM3        1.677613e-195   3.879339 0.975 0.049 2.300678e-191
## CFD           2.401156e-193   3.405492 0.938 0.038 3.292945e-189
## CD68          2.900384e-191   3.020484 0.926 0.035 3.977587e-187
## RP11-290F20.3 2.513244e-186   2.720057 0.840 0.017 3.446663e-182
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
## # A tibble: 18 × 7
## # Groups:   cluster [9]
##        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene    
##        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>   
##  1 9.57e- 88       1.36 0.447 0.108 1.31e- 83 0       CCR7    
##  2 3.75e-112       1.09 0.912 0.592 5.14e-108 0       LDHB    
##  3 0               5.57 0.996 0.215 0         1       S100A9  
##  4 0               5.48 0.975 0.121 0         1       S100A8  
##  5 1.06e- 86       1.27 0.981 0.643 1.45e- 82 2       LTB     
##  6 2.97e- 58       1.23 0.42  0.111 4.07e- 54 2       AQP3    
##  7 0               4.31 0.936 0.041 0         3       CD79A   
##  8 9.48e-271       3.59 0.622 0.022 1.30e-266 3       TCL1A   
##  9 5.61e-202       3.10 0.983 0.234 7.70e-198 4       CCL5    
## 10 7.25e-165       3.00 0.577 0.055 9.95e-161 4       GZMK    
## 11 3.51e-184       3.31 0.975 0.134 4.82e-180 5       FCGR3A  
## 12 2.03e-125       3.09 1     0.315 2.78e-121 5       LST1    
## 13 3.13e-191       5.32 0.961 0.131 4.30e-187 6       GNLY    
## 14 7.95e-269       4.83 0.961 0.068 1.09e-264 6       GZMB    
## 15 1.48e-220       3.87 0.812 0.011 2.03e-216 7       FCER1A  
## 16 1.67e- 21       2.87 1     0.513 2.28e- 17 7       HLA-DPB1
## 17 1.92e-102       8.59 1     0.024 2.63e- 98 8       PPBP    
## 18 9.25e-186       7.29 1     0.011 1.27e-181 8       PF4
##Seurat has several tests for differential expression which can be set 
##with the test.use parameter (see our DE vignette for details). For example, 
##the ROC test returns the ‘classification power’ for any individual marker 
##(ranging from 0 - random, to 1 - perfect).

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

##We include several tools for visualizing marker expression. VlnPlot() (shows expression probability distributions 
##across clusters), and FeaturePlot() (visualizes feature expression on a tSNE or PCA plot) are our most commonly 
##used visualizations. We also suggest exploring RidgePlot(), CellScatter(), and DotPlot() as additional methods to view your dataset.

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))


# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)


FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))


##DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are 
##plotting the top 20 markers (or all markers if less than 20) for each cluster.

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


##Assigning cell type identity to clusters
##Fortunately in the case of this dataset, we can use canonical 
##markers to easily match the unbiased clustering to known cell types:

##  Cluster ID	Markers	Cell Type
##0	IL7R, CCR7	Naive CD4+ T
##1	CD14, LYZ	CD14+ Mono
##2	IL7R, S100A4	Memory CD4+
##  3	MS4A1	B
##4	CD8A	CD8+ T
##5	FCGR3A, MS4A7	FCGR3A+ Mono
##6	GNLY, NKG7	NK
##7	FCER1A, CST3	DC
##8	PPBP	Platelet
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


saveRDS(pbmc, file = "/oak/stanford/groups/cgawad/R_Testing/Seurat/filtered_gene_bc_matrices/hg19/pbmc3k_final.rds")


##############################################
seuratToURD2 <- function(seurat.object) {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    # Create an empty URD object
    ds <- new("URD")
    
    # Copy over data
    ds@logupx.data <- as(as.matrix(seurat.object@assays$RNA@data), "dgCMatrix")
    if(!any(dim(seurat.object@assays$RNA@counts) == 0)) ds@count.data <- as(as.matrix(seurat.object@assays$RNA@counts[rownames(seurat.object@assays$RNA@data), colnames(seurat.object@assays$RNA@data)]), "dgCMatrix")
    
    # Copy over metadata
    ## TO DO - grab kmeans clustering info
    get.data <- NULL
    if (.hasSlot(seurat.object, "data.info")) { 
      get.data <- as.data.frame(seurat.object@assays$RNA@data.info)
    } else if (.hasSlot(seurat.object, "meta.data")) { 
      get.data <- as.data.frame(seurat.object@meta.data) 
    }
    if(!is.null(get.data)) {
      di <- colnames(get.data)
      m <- grep("res|cluster|Res|Cluster", di, value=T, invert = T) # Put as metadata if it's not the result of a clustering.
      discrete <- apply(get.data, 2, function(x) length(unique(x)) / length(x))
      gi <- di[which(discrete <= 0.015)]
      ds@meta <- get.data[,m,drop=F]
      ds@group.ids <- get.data[,gi,drop=F]
    }
    
    # Copy over var.genes
    if(length(seurat.object@assays$RNA@var.features > 0)) ds@var.genes <- seurat.object@assays$RNA@var.features
    
    # Move over tSNE projection
    if (.hasSlot(seurat.object, "tsne.rot")) {
      if(!any(dim(seurat.object@tsne.rot) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@tsne.rot)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("tsne" %in% names(seurat.object@reductions)) && !any(dim(seurat.object@reductions$tsne) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@reductions$tsne@cell.embeddings)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    }
    
    # Move over PCA results
    if (.hasSlot(seurat.object, "pca.x")) {
      if(!any(dim(seurat.object@pca.x) == 0)) {
        ds@pca.load <- seurat.object@pca.x
        ds@pca.scores <- seurat.object@pca.rot
        warning("Need to set which PCs are significant in @pca.sig")
      }
      ## TO DO: Convert SVD to sdev
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("pca" %in% names(seurat.object@reductions)) && !any(dim(Loadings(seurat.object, reduction = "pca")) == 0)) {
        ds@pca.load <- as.data.frame(Loadings(seurat.object, reduction = "pca"))
        ds@pca.scores <- as.data.frame(seurat.object@reductions$pca@cell.embeddings)
        ds@pca.sdev <- seurat.object@reductions$pca@stdev
        ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
      }
    }
    return(ds)
  } else {
    stop("Package Seurat is required for this function. To install: install.packages('Seurat')\n")
  }
}

pbmc_urd<-seuratToURD2(pbmc)

pbmc_urd <- calcPCA(pbmc_urd, mp.factor = 2)

pcSDPlot(pbmc_urd)

set.seed(19)

pbmc_urd <- calcTsne(object = pbmc_urd )

plotDim(pbmc_urd , "seurat_clusters", plot.title = "tSNE: seurat_clusters")

plotDim(pbmc_urd , "CD4", plot.title="tSNE: noto expression (CD4)")


pbmc_urd <- calcDM(pbmc_urd, knn = 100, sigma=16)

plotDimArray(pbmc_urd, reduction.use = "dm", dims.to.plot = 1:8, outer.title = "Diffusion Map (Sigma 16, 100 NNs): Clusters", label="seurat_clusters", plot.title="", legend=F)

plotDim(pbmc_urd, "seurat_clusters", transitions.plot = 10000, plot.title="Developmental stage (with transitions)")



#Here we use all cells from the first stage as the root
root.cells <- cellsInCluster(pbmc_urd, "seurat_clusters", "0")

# Then we run 'flood' simulations
pbmc_urd.floods <- floodPseudotime(pbmc_urd, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)

# The we process the simulations into a pseudotime
pbmc_urd <- floodPseudotimeProcess(pbmc_urd, pbmc_urd.floods, floods.name="pseudotime")

pseudotimePlotStabilityOverall(pbmc_urd)

plotDim(pbmc_urd, "pseudotime")

plotDists(pbmc_urd, "pseudotime", "seurat_clusters", plot.title="Pseudotime by seurat_clusters")


# Create a subsetted object of just those cells from the final stage
axial.6somite <- urdSubset(pbmc_urd, cells.keep=cellsInCluster(pbmc_urd, "seurat_clusters", "5"))

# Get variable genes for each group of 3 stages
# (Normally would do this for each stage, but there are not very many cells in this subset of the data)
# diffCV.cutoff can be varied to include more or fewer genes.
stages <- sort(unique((pbmc_urd@group.ids$seurat_clusters)))
stages              

var.by.stage <- lapply(seq(1,9,1), function(n) {
  findVariableGenes(pbmc_urd, cells.fit=cellsInCluster(pbmc_urd, "seurat_clusters", stages[n]), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, main.use=paste0("Stages ", stages[n]), do.plot=T)
})

var.genes <- sort(unique(unlist(var.by.stage)))
pbmc_urd@var.genes <- var.genes

axial.6somite@var.genes <- var.by.stage[[5]]

# Use the variable genes that were calculated only on the final group of stages (which contain the last stage).
# Calculate PCA and tSNE
axial.6somite <- calcPCA(axial.6somite, mp.factor = 1)
pcSDPlot(axial.6somite)

set.seed(20)
axial.6somite <- calcTsne(axial.6somite)

# Calculate graph clustering of these cells
axial.6somite <- graphClustering(axial.6somite, num.nn = 50, do.jaccard=T, method="Louvain")
plotDim(axial.6somite, "Louvain-50", plot.title = "Louvain (50 NN) graph clustering", point.size=3)
plotDim(axial.6somite, "TNFRSF9", plot.title="HE1A (Differentiated prechordal plate marker)")


pbmc_urd@group.ids[rownames(axial.6somite@group.ids), "tip.clusters"] <- axial.6somite@group.ids$`Louvain-50`
pbmc_urd.ptlogistic <- pseudotimeDetermineLogistic(pbmc_urd, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)
plotDim(pbmc_urd, "tip.clusters", plot.title="Cells in each tip")

plotDim(pbmc_urd, "visitfreq.log.1", plot.title="Visitation frequency from tip 1 (log10)", transitions.plot=10000)




pbmc_urd.tree <- loadTipCells(pbmc_urd, "tip.clusters")   

pbmc_urd.tree <- buildTree(pbmc_urd.tree, pseudotime = "pseudotime", tips.use=1:2,
                           divergence.method = "preference", cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)
