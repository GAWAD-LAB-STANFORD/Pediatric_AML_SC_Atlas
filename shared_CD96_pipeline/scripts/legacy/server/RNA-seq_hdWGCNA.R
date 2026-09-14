# ====================================================================
# RNA-seq_hdWGCNA.R  |  Figures S6/S7 (hdWGCNA co-expression modules)
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : ALSF_RNA-seq_for_CytoTrace_HBM_FAB_phenotype.csv, ALSF_RNA-seq_for_CytoTrace_HBM_expr.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript RNA-seq_hdWGCNA.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

library(dplyr)
library(Seurat)
library(patchwork)
library(devtools)
library(URD)
use_condaenv(condaenv = "scRNA", conda = "$OAK/Sequencing_Analysis_Tools/scRNA-Python/etc/profile.d/conda.sh")


install.packages("reticulate")
BiocManager::install("sva")
devtools::install_local("CytoTRACE_0.3.3.tar.gz")
BiocManager::install("WGCNA") 
devtools::install_github('smorabit/hdWGCNA', ref='dev')
library(reticulate)
library(CytoTRACE)
#################read csv data to R dataframe
ALSF_HBM2 <-read.delim2("/oak/stanford/groups/cgawad/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/Seurat/ALSF_RNA-seq_for_CytoTrace_HBM_expr.csv",sep=',')
ALSF_HBM2_pheno<-read.delim2("/oak/stanford/groups/cgawad/Cancer_Studies/SC_RNA_SEQ/ALSF_AML/scanpy/Seurat/ALSF_RNA-seq_for_CytoTrace_HBM_FAB_phenotype.csv",sep=',')
ALSF_HBM2__mtx<-ALSF_HBM2[,-1]
rownames(ALSF_HBM2__mtx)<-ALSF_HBM2[,1]
ALSF_HBM2__mtx[,1:ncol(ALSF_HBM2__mtx)]<-lapply(ALSF_HBM2__mtx[,1:ncol(ALSF_HBM2__mtx)], as.numeric)
################run cytoTrace
results <- CytoTRACE(ALSF_HBM2__mtx, subsamplesize = 1000)
###################plot
plotCytoTRACE(results, phenotype = marrow_10x_pheno, gene = "CD34")
