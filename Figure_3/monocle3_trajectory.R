#!/usr/bin/env Rscript
# Monocle3 trajectory of the 49 leukemic states + 8 normal compartments (subsampled). Batch-aligned
# by patient (SampleID), single principal graph, ordered from a normal-HSPC root. Outputs a UMAP with
# the trajectory coloured by Monocle pseudotime and per-cell/per-group pseudotime for ordering the
# states with normal anchors. Independent cross-check of the diffusion-pseudotime ordering. No fabrication.
suppressPackageStartupMessages({library(monocle3); library(Matrix); library(dplyr); library(ggplot2); library(igraph)})
IN  <- "/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad/monocle_input"
OUT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/04_Main_Figures/Figure_3_panels"
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
set.seed(1)
counts <- as(Matrix::readMM(file.path(IN,"counts.mtx")), "CsparseMatrix")   # genes x cells
cells  <- read.csv(file.path(IN,"cells.csv")); genes <- read.csv(file.path(IN,"genes.csv"))
rownames(counts) <- make.unique(genes$gene); colnames(counts) <- cells$barcode
gm <- data.frame(gene_short_name=rownames(counts), row.names=rownames(counts))
cm <- data.frame(cells, row.names=cells$barcode)
cds <- new_cell_data_set(counts, cell_metadata=cm, gene_metadata=gm)
cds <- preprocess_cds(cds, num_dim=50)
cds <- align_cds(cds, alignment_group="SampleID")     # correct patient/batch effects
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition=FALSE)           # single graph across all cells
# root = principal-graph node with the most normal-HSPC cells
cv <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
cv <- as.matrix(cv[colnames(cds), , drop=FALSE])
hspc <- which(colData(cds)$group == "NORM_HSPC")
root <- V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(cv[hspc,1]))))]
cds <- order_cells(cds, root_pr_nodes=root)
pt <- pseudotime(cds); pt[!is.finite(pt)] <- NA
colData(cds)$monocle_pt <- pt
out <- data.frame(barcode=colnames(cds), group=colData(cds)$group, monocle_pt=pt, dpt=colData(cds)$dpt)
write.csv(out, file.path(D3,"monocle_pseudotime.csv"), row.names=FALSE)
# export Monocle UMAP coordinates (for the custom LS19-highlight panel)
umap <- SingleCellExperiment::reducedDims(cds)[["UMAP"]]
coords <- data.frame(barcode=colnames(cds), UMAP1=umap[,1], UMAP2=umap[,2],
                     group=colData(cds)$group, monocle_pt=pt)
write.csv(coords, file.path(D3,"monocle_umap_coords.csv"), row.names=FALSE)
# export the learned principal graph (UMAP frame) + per-node pseudotime, for drawing the
# directed trajectory (arrows oriented low->high pseudotime) on the custom LS19-highlight panel
pg <- principal_graph(cds)[["UMAP"]]
nm <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
gnodes <- data.frame(node=rownames(nm), UMAP1=nm[,1], UMAP2=nm[,2])
cell_node <- igraph::V(pg)$name[cv[,1]]
node_pt <- tapply(pt, cell_node, function(z) median(z, na.rm=TRUE))
gnodes$pt <- as.numeric(node_pt[gnodes$node])
el <- igraph::as_edgelist(pg)
write.csv(gnodes, file.path(D3,"monocle_graph_nodes.csv"), row.names=FALSE)
write.csv(data.frame(from=el[,1], to=el[,2]), file.path(D3,"monocle_graph_edges.csv"), row.names=FALSE)
cat("exported principal graph:", nrow(gnodes), "nodes,", nrow(el), "edges; root =", root, "\n")
cat("root node:", root, "| pseudotime range:", round(range(pt,na.rm=TRUE),3), "| NA(unreachable):", sum(is.na(pt)), "\n")
# UMAP with trajectory coloured by pseudotime
p1 <- plot_cells(cds, color_cells_by="pseudotime", label_cell_groups=FALSE, label_leaves=FALSE,
                 label_branch_points=FALSE, label_roots=FALSE, cell_size=0.5, trajectory_graph_color="grey30")
ggsave(file.path(OUT,"Figure_3_monocle_umap.png"), p1, width=6.6, height=5.6, dpi=170)
ggsave(file.path(OUT,"Figure_3_monocle_umap.pdf"), p1, width=6.6, height=5.6)
# normal-compartment sanity + save medians
med <- out %>% filter(!is.na(monocle_pt)) %>% group_by(group) %>% summarise(med=median(monocle_pt), n=n(), .groups="drop")
cat("normal compartment median Monocle pseudotime (HSPC should be lowest):\n")
print(med %>% filter(grepl("^NORM_",group)) %>% arrange(med) %>% as.data.frame())
cat("wrote monocle_pseudotime.csv and Figure_3_monocle_umap\n")
