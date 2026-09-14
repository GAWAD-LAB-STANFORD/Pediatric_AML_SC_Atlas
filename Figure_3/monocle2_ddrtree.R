#!/usr/bin/env Rscript
# Monocle2 DDRTree trajectory of the 49 leukemic states + 8 normal compartments. DDRTree is O(n)-heavy,
# so cells are subsampled (<=60/group) to keep the reverse-graph embedding tractable. Patient/batch is
# regressed out in reduceDimension (residualModelFormulaStr = ~SampleID), and the tree is rooted at the
# State containing the most normal HSPC. Exports the DDRTree cell coordinates + tree skeleton for the
# custom LS19-highlight panel. This is the classic branched-"line" trajectory. No fabricated values.
suppressPackageStartupMessages({library(monocle); library(Matrix); library(igraph)})
# monocle 2.22.0's orderCells() -> extract_ddrtree_ordering() calls igraph::graph.dfs(..., neimode=),
# which modern igraph made defunct. Rewrite just that one defunct call to the current dfs(mode=) API,
# returning order/father as numeric vertex indices (the contract graph.dfs provided). No igraph downgrade.
.eddo <- monocle:::extract_ddrtree_ordering
.txt  <- paste(deparse(.eddo), collapse = "\n")
.txt  <- gsub("graph\\.dfs\\s*\\(\\s*dp_mst\\s*,\\s*root\\s*=\\s*root_cell\\s*,\\s*neimode\\s*=\\s*\"all\"\\s*,\\s*unreachable\\s*=\\s*FALSE\\s*,\\s*father\\s*=\\s*TRUE\\s*\\)",
              "local({ .r <- igraph::dfs(dp_mst, root = root_cell, mode = \"all\", unreachable = FALSE, father = TRUE); list(order = as.integer(.r$order), father = as.integer(.r$father)) })",
              .txt)
.eddo2 <- eval(parse(text = .txt)); environment(.eddo2) <- asNamespace("monocle")
assignInNamespace("extract_ddrtree_ordering", .eddo2, ns = "monocle")
if (grepl("graph.dfs", paste(deparse(monocle:::extract_ddrtree_ordering), collapse=" ")))
  stop("patch failed: graph.dfs still present in extract_ddrtree_ordering")
cat("patched extract_ddrtree_ordering (graph.dfs -> igraph::dfs mode=)\n")
# project2MST also uses igraph calls defunct in igraph 2.x: nei() -> .nei(), and the deprecated
# graph.adjacency() -> graph_from_adjacency_matrix(). Rewrite both.
.p2m  <- paste(deparse(monocle:::project2MST), collapse = "\n")
.p2m  <- gsub("\\bnei\\(", ".nei(", .p2m)
.p2m  <- gsub("graph\\.adjacency\\(", "igraph::graph_from_adjacency_matrix(", .p2m)
.p2m  <- gsub("class\\(projection\\)\\s*!=\\s*\"matrix\"", "!is.matrix(projection)", .p2m)  # R>=4.2 class() length
.p2mf <- eval(parse(text = .p2m)); environment(.p2mf) <- asNamespace("monocle")
assignInNamespace("project2MST", .p2mf, ns = "monocle")
cat("patched project2MST (nei -> .nei, graph.adjacency -> graph_from_adjacency_matrix)\n")
IN  <- "/private/tmp/claude-503/-Users-chuckgawad-Desktop-ALSF-AML-2026-updated/f5e93d65-31d5-411e-b754-6f7c504928b7/scratchpad/monocle_input"
D3  <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/__SUBMISSION_PACKAGE_v2_recluster/08_Source_Data/Figure_3/v3"
set.seed(1)
counts <- as(Matrix::readMM(file.path(IN,"counts.mtx")), "CsparseMatrix")   # genes x cells
cells  <- read.csv(file.path(IN,"cells.csv")); genes <- read.csv(file.path(IN,"genes.csv"))
rownames(counts) <- make.unique(genes$gene); colnames(counts) <- cells$barcode
# subsample <=60 cells/group for DDRTree tractability
idx <- unlist(lapply(split(seq_len(nrow(cells)), cells$group),
                     function(ii) if (length(ii) <= 60) ii else sample(ii, 60)))
counts <- counts[, idx]; cells <- cells[idx, ]
rownames(cells) <- cells$barcode
cat("subsampled to", ncol(counts), "cells across", length(unique(cells$group)), "groups\n")
pd <- new("AnnotatedDataFrame", data = cells)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(counts),
                                                  row.names = rownames(counts)))
cds <- newCellDataSet(counts, phenoData = pd, featureData = fd,
                      lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
# NOTE: monocle 2.22.0's estimateDispersions() calls the now-defunct dplyr::group_by_(); rather than
# downgrade dplyr (which the rest of the environment needs), select ordering genes manually from the
# size-factor-normalised counts (top overdispersed genes). DDRTree with norm_method="log" does not need
# the dispersion model.
X  <- exprs(cds); sf <- sizeFactors(cds)
Xn <- X %*% Matrix::Diagonal(x = 1 / sf)                 # normalise columns (sparse)
gm  <- Matrix::rowMeans(Xn); gm2 <- Matrix::rowMeans(Xn * Xn)
gv  <- (gm2 - gm^2) * (ncol(Xn) / (ncol(Xn) - 1))
disp <- gv / (gm^2)
ok  <- gm >= 0.1 & is.finite(disp)
ordering <- names(sort(disp[ok], decreasing = TRUE))[seq_len(min(1500, sum(ok)))]
cat("ordering genes:", length(ordering), "\n")
cds <- setOrderingFilter(cds, ordering)
cds <- reduceDimension(cds, max_components = 2, reduction_method = "DDRTree", norm_method = "log",
                       residualModelFormulaStr = "~SampleID")
cds <- orderCells(cds)
# root at the State with the most normal-HSPC cells
hspc_by_state <- tapply(pData(cds)$group == "NORM_HSPC", pData(cds)$State, sum)
root_state <- names(which.max(hspc_by_state))
cds <- orderCells(cds, root_state = root_state)
cat("root State:", root_state, "| pseudotime range:", round(range(pData(cds)$Pseudotime), 2), "\n")
# export cell coordinates
S <- t(reducedDimS(cds))
cellout <- data.frame(barcode = colnames(cds), DDR1 = S[,1], DDR2 = S[,2],
                      Pseudotime = pData(cds)$Pseudotime, State = pData(cds)$State,
                      group = pData(cds)$group)
write.csv(cellout, file.path(D3, "ddrtree_cells.csv"), row.names = FALSE)
# export tree skeleton (principal graph node coords + edges)
K <- t(reducedDimK(cds))
nodes <- data.frame(node = rownames(K), K1 = K[,1], K2 = K[,2])
mst <- minSpanningTree(cds)
el <- igraph::as_edgelist(mst)
write.csv(nodes, file.path(D3, "ddrtree_nodes.csv"), row.names = FALSE)
write.csv(data.frame(from = el[,1], to = el[,2]), file.path(D3, "ddrtree_edges.csv"), row.names = FALSE)
# diagnostics: median pseudotime by compartment (HSPC should be lowest)
med <- tapply(pData(cds)$Pseudotime, pData(cds)$group, median)
cat("median pseudotime (HSPC lowest expected), key groups:\n")
print(round(sort(med)[c("NORM_HSPC","LS_19","NORM_DC","NORM_Monocyte","NORM_T_NK")], 2))
cat("LS19 median pseudotime:", round(median(pData(cds)$Pseudotime[pData(cds)$group=="LS_19"]), 2),
    "| overall max:", round(max(pData(cds)$Pseudotime), 2), "\n")
cat("wrote ddrtree_cells.csv / ddrtree_nodes.csv / ddrtree_edges.csv\n")
