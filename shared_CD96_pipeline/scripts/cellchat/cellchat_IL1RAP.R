# ====================================================================
# cellchat_IL1RAP.R  |  CD96 figure pipeline component
# reproducibility header (auto-generated; edit freely below this block)
# Inputs : cellchat_IL1_input.csv
# Outputs: (writes to ../figures or ../source_data — see body)
# Run    : Rscript cellchat_IL1RAP.R   (or: bash run_all.sh)  |  Deps: see REPRODUCIBILITY_CHECKLIST.md
# Source data for this figure is in __SUBMISSION_PACKAGE/08_Source_Data/.
# ====================================================================

# CellChat IL-1 pathway (IL1A/IL1B -> IL1R1_IL1RAP; IL36A/IL36G -> IL1RL2_IL1RAP) in
# pediatric AML — the IL1RAP axis. IL1RAP is a RECEPTOR co-subunit, so AML is the
# RECEIVER (mirror-image of the CD96 panel, where AML is the sender). This captures the
# microenvironment-IL1 -> AML-LSC maintenance loop. Inferred WITHIN each sample; aggregated
# over the IL1RAP-targetable samples (>=10% IL1RAP+ blasts). truncatedMean, trim=0.
suppressPackageStartupMessages({ library(CellChat); library(Matrix) })
HERE <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/AML_2026_Alternate/scripts/cellchat"
setwd(HERE)
ROOT <- "/Users/chuckgawad/Desktop/ALSF_AML_2026_updated/AML_2026_Alternate"
SD <- file.path(ROOT, "source_data"); ASSETS <- file.path(ROOT, "assets")
MINC <- 10; POS_MIN <- 10                  # >=10% IL1RAP+ AML cells = targetable (receiver)
GENES <- c("IL1A", "IL1B", "IL36A", "IL36G", "IL1R1", "IL1RL1", "IL1RL2", "IL1RAP")
d <- read.csv("cellchat_IL1_input.csv", stringsAsFactors = FALSE, check.names = FALSE)

aml <- d[d$compartment == "AML (leukemic)", ]
frac <- tapply(aml$IL1RAP > 0, aml$sample, mean) * 100
sampinfo <- data.frame(sample = names(frac), pct_il1rap_blasts = as.numeric(frac))
sampinfo$targetable <- sampinfo$pct_il1rap_blasts >= POS_MIN
targetable <- sampinfo$sample[sampinfo$targetable]
cat("IL1RAP-targetable samples:", length(targetable), "of", nrow(sampinfo), "\n")
write.csv(sampinfo, file.path(SD, "il1rap_sample_info.csv"), row.names = FALSE)

DBuse <- subsetDB(CellChatDB.human, search = "IL1", key = "pathway_name")
# keep IL1RAP-containing receptor complexes only (the IL1R2 decoy doesn't signal and
# its gene isn't in the matrix), with ligands present in the atlas
DBuse$interaction <- DBuse$interaction[DBuse$interaction$ligand %in% GENES &
                                       grepl("IL1RAP", DBuse$interaction$receptor), ]
LRtab <- DBuse$interaction[, c("interaction_name", "ligand", "receptor")]
allcomp <- sort(unique(d$compartment)); LRn <- LRtab$interaction_name
z <- function() array(0, c(length(allcomp), length(allcomp), length(LRn)),
                      dimnames = list(allcomp, allcomp, LRn))
probsum <- z(); assess <- z(); sigcnt <- z()
persample <- list()

for (s in sampinfo$sample) {
  ds <- d[d$sample == s, ]
  keep <- names(which(table(ds$compartment) >= MINC))
  if (length(keep) < 2) next
  ds <- ds[ds$compartment %in% keep, ]
  mat <- t(as.matrix(ds[, GENES])); colnames(mat) <- ds$cell; rownames(mat) <- GENES
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
  }, error = function(e) { cat("  [skip]", s, "::", conditionMessage(e), "\n"); NULL })
  if (is.null(cc)) next
  pr <- cc@net$prob; pv <- cc@net$pval; g <- levels(cc@idents)
  istgt <- s %in% targetable
  for (L in dimnames(pr)[[3]]) for (a in g) for (b in g) {
    persample[[length(persample) + 1]] <- data.frame(
      sample = s, targetable = istgt, source = a, target = b, interaction = L,
      prob = pr[a, b, L], pval = pv[a, b, L])
    if (istgt) {
      assess[a, b, L]  <- assess[a, b, L] + 1
      probsum[a, b, L] <- probsum[a, b, L] + pr[a, b, L]
      if (!is.na(pv[a, b, L]) && pv[a, b, L] <= 0.05) sigcnt[a, b, L] <- sigcnt[a, b, L] + 1
    }
  }
}
ps <- do.call(rbind, persample)
write.csv(ps, file.path(SD, "il1rap_persample_long.csv"), row.names = FALSE)
meanprob <- probsum / assess; meanprob[!is.finite(meanprob)] <- 0
fracsig  <- sigcnt / assess;  fracsig[!is.finite(fracsig)] <- 0

rows <- list()
for (L in LRn) for (a in allcomp) for (b in allcomp) if (assess[a, b, L] > 0)
  rows[[length(rows) + 1]] <- data.frame(source = a, target = b, interaction = L,
    ligand = LRtab$ligand[LRtab$interaction_name == L],
    receptor = LRtab$receptor[LRtab$interaction_name == L],
    mean_prob = meanprob[a, b, L], n_assessable = assess[a, b, L],
    n_sig = sigcnt[a, b, L], frac_sig = fracsig[a, b, L])
cat(sprintf("\nDIAG: persample rows=%d, targetable rows=%d, sum(assess)=%d, aggregated rows=%d\n",
            nrow(ps), sum(ps$targetable), sum(assess), length(rows)))
agg <- if (length(rows)) do.call(rbind, rows) else data.frame()
if (nrow(agg)) {
  agg <- agg[order(-agg$mean_prob), ]
  write.csv(agg, file.path(SD, "il1rap_communication_aggregated.csv"), row.names = FALSE)
  cat("=== top interactions signalling TO AML (IL1RAP on leukemic cells) ===\n")
  toaml <- agg[agg$target == "AML (leukemic)",
               c("source", "interaction", "mean_prob", "n_assessable", "n_sig", "frac_sig")]
  print(head(toaml, 12), row.names = FALSE)
} else {
  cat("WARNING: no aggregated interactions (assess all zero).\n")
}

## circle: senders -> AML (AML as receiver), node size = #cells in targetable samples
dt <- d[d$sample %in% targetable, ]
vsz <- as.numeric(table(factor(dt$compartment, levels = allcomp)))
W <- apply(meanprob, c(1, 2), sum)
Wrecv <- W; Wrecv[, colnames(Wrecv) != "AML (leukemic)"] <- 0   # edges INTO AML only
png(file.path(ASSETS, "fig6_il1rap_circle.png"), 6.4, 6.4, units = "in", res = 220, bg = "white")
netVisual_circle(Wrecv, vertex.weight = vsz, weight.scale = TRUE, label.edge = FALSE,
                 title.name = "IL-1 Signalling TO AML Blasts (IL1RAP+ Samples)")
dev.off()
cat("\nDone -> il1rap_{aggregated,persample_long,sample_info}.csv + assets/fig6_il1rap_circle.png\n")
