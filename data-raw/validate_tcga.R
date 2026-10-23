## Validate thyroidBRS against the published TCGA-THCA BRAF-RAS Score.
##
## The point of this script is to measure agreement with something that did
## NOT go into the fit. The centroids are fitted on driver-mutation status
## (BRAF-V600E vs RAS); the published `BRAF_RAF_score` / `BRAF_RAF_class`
## are only ever used to check the result. Fitting on the published class
## instead is circular — it is the sign of the published score.
##
##   THCA_LOG2TPM=<thca_log2tpm.rds> Rscript data-raw/validate_tcga.R
##
## THCA_LOG2TPM is the log2(TPM + 1) matrix that data-raw/build_inputs.R
## writes: TCGA-THCA primary tumors, gene symbols as row names, patient
## barcodes (TCGA-XX-XXXX) as column names. Build it with that script rather
## than by hand: how duplicated symbols are collapsed and which aliquot is
## kept per patient both move the figures below, and build_inputs.R is where
## those decisions are written down.

## Assisted-by: Claude Opus 5, Claude Fable 5.1 and Claude Opus 5.5
## (Anthropic). See the Provenance section of README.md.

library(thyroidBRS)

source(file.path("data-raw", "record_session.R"))

expr_path <- Sys.getenv("THCA_LOG2TPM", "thca_log2tpm.rds")
if (!file.exists(expr_path)) {
    stop("Set THCA_LOG2TPM to a log2(TPM + 1) matrix; see the header.")
}
expr <- readRDS(expr_path)

ref <- read.csv(
    system.file("extdata", "thca_reference.csv", package = "thyroidBRS"),
    stringsAsFactors = FALSE
)
rownames(ref) <- ref$patient

## ---- fit on mutation status -------------------------------------------
mut <- setNames(ref$driver_group, ref$patient)
mut <- mut[mut %in% c("BRAF_V600E", "RAS") & names(mut) %in% colnames(expr)]

fit <- brs_fit(expr, mut)
print(fit)

preds <- predict(fit, expr)
rownames(preds) <- preds$sample

## ---- validate against the published score -----------------------------
pub <- ref[ref$patient %in% preds$sample & !is.na(ref$published_brs), ]
i <- pub$patient
held_out <- !(i %in% fit$reference_samples)

cat("\nAgreement with the published BRS (n = ", nrow(pub), ")\n", sep = "")
cat("  Spearman vs BRAF_RAF_score   : ",
    round(cor(preds[i, "brs_score"], pub$published_brs,
              method = "spearman"), 4), "\n", sep = "")
## brs_scaled is rescaled across the set being scored, and the published
## values were rescaled across the 391 that had one. The rescaling is linear
## within each sign but not across zero, so this comparison only holds while
## both extremes fall inside the shared subset. They do here; check it before
## carrying the number to another cohort.
cat("  Pearson, rescaled            : ",
    round(cor(preds[i, "brs_scaled"], pub$published_brs), 4), "\n", sep = "")
cat("  Class concordance, all       : ",
    round(100 * mean(preds[i, "brs_class"] == pub$published_class), 1),
    "%\n", sep = "")
cat("  Class concordance, held out  : ",
    round(100 * mean(preds[i, "brs_class"][held_out] ==
                     pub$published_class[held_out]), 1),
    "% (n = ", sum(held_out), ")\n", sep = "")
print(table(published = pub$published_class, predicted = preds[i, "brs_class"]))

## ---- cross-validate on the reference set ------------------------------
set.seed(1)
fold <- sample(rep(seq_len(10), length.out = length(mut)))
cv_fold <- function(k) {
    f <- brs_fit(expr, mut[fold != k])
    in_fold <- names(mut)[fold == k]
    setNames(predict(f, expr[, in_fold, drop = FALSE])$brs_class, in_fold)
}
cv <- unlist(lapply(seq_len(10), cv_fold))[names(mut)]
truth <- ifelse(mut == "BRAF_V600E", "Braf-like", "Ras-like")
cat("\n10-fold CV against mutation status (n = ", length(mut), "): ",
    round(100 * mean(cv == truth), 1), "%\n", sep = "")

## ---- coverage gained ---------------------------------------------------
new <- setdiff(preds$sample, pub$patient)
cat("\nScored ", nrow(preds), " patients; ", length(new),
    " have no published BRS.\n", sep = "")
print(table(preds[new, "brs_class"]))

## ---- orientation check -------------------------------------------------
## Block 1 of Figure S7A should be higher in RAS, block 2 higher in BRAF.
delta <- fit$centroid_ras - fit$centroid_braf
key <- c(setNames(brs_genes$up_in, brs_genes$current_symbol),
         setNames(brs_genes$up_in, brs_genes$original_symbol))
expected <- key[names(delta)]
cat("\nOrientation matches Figure S7A for ",
    sum(ifelse(delta > 0, "RAS", "BRAF") == expected), " of ", length(delta),
    " genes\n", sep = "")

## ---- FLJ23867 as QSOX1 -------------------------------------------------
## FLJ23867 is a clone name whose NCBI Gene record was merged into QSOX1 in
## 2018, so in 2014 it was a separate Gene ID. Check that QSOX1 behaves like
## the rest of its block: BRAF-RAS difference on the reference groups, and
## Spearman correlation with the axis built from the other signature genes
## (block 2 minus block 1 of per-gene z-scores, the gene itself left out).
sig <- fit$genes_used
blk <- setNames(brs_genes$block, brs_genes$current_symbol)
blk <- c(blk, setNames(brs_genes$block, brs_genes$original_symbol))[sig]
z <- t(scale(t(expr[sig, ])))
b1 <- sig[blk == 1L]
b2 <- sig[blk == 2L]
rho_b2 <- vapply(b2, function(g) {
    axis <- colMeans(z[setdiff(b2, g), ]) - colMeans(z[b1, ])
    cor(z[g, ], axis, method = "spearman")
}, numeric(1))
lfc_b2 <- vapply(b2, function(g) {
    mean(expr[g, names(mut)[mut == "BRAF_V600E"]]) -
        mean(expr[g, names(mut)[mut == "RAS"]])
}, numeric(1))
others <- setdiff(b2, "QSOX1")
cat("\nQSOX1 (published as FLJ23867) against the other ", length(others),
    " block-2 genes\n", sep = "")
cat("  log2 BRAF - RAS : ", round(lfc_b2[["QSOX1"]], 2), " (others ",
    paste(round(range(lfc_b2[others]), 2), collapse = " to "), ", median ",
    round(median(lfc_b2[others]), 2), ")\n", sep = "")
cat("  Spearman vs axis: ", round(rho_b2[["QSOX1"]], 3), " (others ",
    paste(round(range(rho_b2[others]), 2), collapse = " to "), ", median ",
    round(median(rho_b2[others]), 3), ")\n", sep = "")

record_session("validate_tcga")
