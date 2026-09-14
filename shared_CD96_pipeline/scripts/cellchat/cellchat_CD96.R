# ====================================================================
# cellchat_CD96.R  |  CD96 figure pipeline component
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : cellchat_input.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript cellchat_CD96.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# CellChat CD96 pathway (CD96->PVR[CD155], CD96->NECTIN1; cell-cell contact) in
# pediatric AML. Focus: CD96 ON THE AML (leukemic) cells = the therapeutic target,
# so AML is the sender (CellChatDB encodes CD96 as the ligand). Communication is
# inferred WITHIN each sample. For the MAIN figure (7E/7F) we aggregate only the
# CD96-TARGETABLE samples (>=10% CD96+ blasts). We ALSO run every sample (incl.
# CD96-cold) and keep per-sample probabilities so the supplement (Fig S8) can show
# the full sender network, per-patient heterogeneity, and the targetable-vs-cold
# contrast. Sparse markers -> per-group mean (truncatedMean, trim=0).
suppressPackageStartupMessages({ library(CellChat); library(Matrix) })
HERE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/cellchat_CD96"; setwd(HERE)
ASSETS <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/CD96_paper_figures_R/assets"
MINC <- 10; CD96POS_MIN <- 10            # >=10% CD96+ AML cells = targetable
d <- read.csv("cellchat_input.csv", stringsAsFactors = FALSE)

## ---- per-sample CD96+ blast fraction; targetable = >=10% CD96+ leukemic ----
aml <- d[d$compartment == "AML (leukemic)", ]
frac <- tapply(aml$CD96 > 0, aml$sample, mean) * 100
nbl  <- tapply(aml$CD96, aml$sample, length)
sampinfo <- data.frame(sample = names(frac), pct_cd96_blasts = as.numeric(frac),
                       n_leukemic = as.numeric(nbl[names(frac)]))
sampinfo$targetable <- sampinfo$pct_cd96_blasts >= CD96POS_MIN
targetable <- sampinfo$sample[sampinfo$targetable]
cat("CD96-targetable samples:", length(targetable), "of", nrow(sampinfo), "\n")
write.csv(sampinfo, "cd96_sample_info.csv", row.names = FALSE)

## ---- autocrine substrate: co-expression on leukemic cells (targetable) ----
acc <- do.call(rbind, lapply(targetable, function(s) {
  a <- aml[aml$sample == s, ]
  data.frame(sample = s, n_leukemic = nrow(a),
             pct_CD96 = mean(a$CD96 > 0) * 100, pct_PVR = mean(a$PVR > 0) * 100,
             pct_NECTIN1 = mean(a$NECTIN1 > 0) * 100,
             pct_CD96_PVR = mean(a$CD96 > 0 & a$PVR > 0) * 100,
             pct_CD96_NECTIN1 = mean(a$CD96 > 0 & a$NECTIN1 > 0) * 100)
}))
write.csv(acc, "cd96_autocrine_coexpr.csv", row.names = FALSE)

DBuse <- subsetDB(CellChatDB.human, search = "CD96", key = "pathway_name")
LRtab <- DBuse$interaction[, c("interaction_name", "ligand", "receptor")]
allcomp <- sort(unique(d$compartment)); LRn <- LRtab$interaction_name
z <- function() array(0, c(length(allcomp), length(allcomp), length(LRn)),
                      dimnames = list(allcomp, allcomp, LRn))
probsum <- z(); assess <- z(); sigcnt <- z()      # aggregated over TARGETABLE only
persample <- list()                               # long per-sample (ALL samples)

for (s in sampinfo$sample) {
  ds <- d[d$sample == s, ]
  keep <- names(which(table(ds$compartment) >= MINC))
  if (length(keep) < 2) next
  ds <- ds[ds$compartment %in% keep, ]
  mat <- t(as.matrix(ds[, c("CD96", "PVR", "NECTIN1")]))
  colnames(mat) <- ds$cell; rownames(mat) <- c("CD96", "PVR", "NECTIN1")
  mat <- as(mat, "CsparseMatrix")
  meta <- data.frame(compartment = factor(ds$compartment, levels = keep), row.names = ds$cell)
  cc <- tryCatch({
    x <- createCellChat(object = mat, meta = meta, group.by = "compartment")
    x@DB <- DBuse; x <- subsetData(x)
    x <- tryCatch(identifyOverExpressedGenes(x),
                  error = function(e) { x@var.features$features <- rownames(mat); x })
    x <- tryCatch(identifyOverExpressedInteractions(x), error = function(e) x)
    x@LR$LRsig <- DBuse$interaction
    x <- computeCommunProb(x, type = "truncatedMean", trim = 0, raw.use = TRUE,
                           population.size = TRUE)
    filterCommunication(x, min.cells = MINC)
  }, error = function(e) NULL)
  if (is.null(cc)) next
  pr <- cc@net$prob; pv <- cc@net$pval; g <- levels(cc@idents)
  istgt <- s %in% targetable
  for (L in dimnames(pr)[[3]]) for (a in g) for (b in g) {
    persample[[length(persample) + 1]] <- data.frame(
      sample = s, targetable = istgt, source = a, target = b, interaction = L,
      prob = pr[a, b, L], pval = pv[a, b, L])
    if (istgt) {                                   # aggregate targetable for 7E/7F
      assess[a, b, L]  <- assess[a, b, L] + 1
      probsum[a, b, L] <- probsum[a, b, L] + pr[a, b, L]
      if (!is.na(pv[a, b, L]) && pv[a, b, L] <= 0.05) sigcnt[a, b, L] <- sigcnt[a, b, L] + 1
    }
  }
}
ps <- do.call(rbind, persample)
write.csv(ps, "cd96_persample_long.csv", row.names = FALSE)

meanprob <- probsum / assess; meanprob[!is.finite(meanprob)] <- 0
fracsig  <- sigcnt / assess;  fracsig[!is.finite(fracsig)] <- 0

## ---- tidy aggregated table (TARGETABLE) -> Fig 7E/7F ----
rows <- list()
for (L in LRn) for (a in allcomp) for (b in allcomp) if (assess[a, b, L] > 0)
  rows[[length(rows) + 1]] <- data.frame(source = a, target = b, interaction = L,
    ligand = LRtab$ligand[LRtab$interaction_name == L],
    receptor = LRtab$receptor[LRtab$interaction_name == L],
    mean_prob = meanprob[a, b, L], n_assessable = assess[a, b, L],
    n_sig = sigcnt[a, b, L], frac_sig = fracsig[a, b, L])
agg <- do.call(rbind, rows); agg <- agg[order(-agg$mean_prob), ]
write.csv(agg, "cd96_communication_aggregated.csv", row.names = FALSE)
cat("=== top AML-source CD96 interactions (CD96 on leukemic cells) ===\n")
print(head(agg[agg$source == "AML (leukemic)",
               c("target","interaction","mean_prob","n_assessable","n_sig","frac_sig")], 12),
      row.names = FALSE)

## ---- circles (node size = #cells in targetable samples) ----
dt  <- d[d$sample %in% targetable, ]
vsz <- as.numeric(table(factor(dt$compartment, levels = allcomp)))
W <- apply(meanprob, c(1, 2), sum)               # sum the 2 L-R pairs
mkcircle <- function(M, file, title) {
  png(file, 6.4, 6.4, units = "in", res = 220, bg = "white")
  netVisual_circle(M, vertex.weight = vsz, weight.scale = TRUE, label.edge = FALSE,
                   title.name = title)
  dev.off()
}
Waml <- W; Waml[rownames(Waml) != "AML (leukemic)", ] <- 0   # AML as sender only
mkcircle(Waml, file.path(ASSETS, "fig7E_cellchat_circle.png"),
         "CD96 Signalling From AML Blasts (CD96+ Samples)")
mkcircle(W, file.path(ASSETS, "figS8A_cellchat_circle_all.png"),
         "CD96 Communication — All Senders (CD96+ Samples)")
mkcircle(W, "cd96_circle_full.png", "CD96 Communication — All Senders (CD96+ Samples)")
cat("\nDone -> aggregated + persample_long + autocrine_coexpr + sample_info; circles in assets/\n")
