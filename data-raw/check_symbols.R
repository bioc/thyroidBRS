## Re-check the BRS signature symbols against current annotation.
## Run after an org.Hs.eg.db update; it prints what, if anything, has drifted.
##
##   Rscript data-raw/check_symbols.R

## Assisted-by: Claude Opus 5, Claude Fable 5.1 and Claude Opus 5.5
## (Anthropic). See the Provenance section of README.md.

library(org.Hs.eg.db)

source(file.path("data-raw", "record_session.R"))
devtools::load_all(".", quiet = TRUE)

approved <- AnnotationDbi::keys(org.Hs.eg.db, "SYMBOL")

resolve <- function(symbol) {
    hit <- tryCatch(
        suppressMessages(
            AnnotationDbi::select(org.Hs.eg.db, keys = symbol,
                                  keytype = "ALIAS", columns = "SYMBOL")
        ),
        error = function(e) NULL
    )
    if (is.null(hit)) return(NA_character_)
    paste(unique(hit$SYMBOL), collapse = "/")
}

cat("org.Hs.eg.db: ", as.character(packageVersion("org.Hs.eg.db")), "\n",
    "checked on:   ", format(Sys.Date()), "\n\n", sep = "")

## 1. Every current_symbol we ship should be an approved symbol today.
stale <- with(brs_genes, current_symbol[!is.na(current_symbol) &
                                        !current_symbol %in% approved])
if (length(stale)) {
    cat("current_symbol no longer approved:\n")
    cat(sprintf("  %-12s -> %s\n", stale,
                vapply(stale, resolve, character(1))), sep = "")
} else {
    cat("All current_symbol values are approved HGNC symbols.\n")
}

## 2. Every original_symbol that is stale should have a mapping.
cat("\nOriginal symbols that are not approved symbols today:\n")
old <- with(brs_genes, original_symbol[!original_symbol %in% approved])
shipped <- brs_genes$current_symbol[match(old, brs_genes$original_symbol)]
cat(sprintf("  %-12s org.Hs.eg.db: %-10s  shipped: %s\n",
            old, vapply(old, resolve, character(1)),
            ifelse(is.na(shipped), "<NA>", shipped)), sep = "")

## 3. Clone names. org.Hs.eg.db carries aliases but not withdrawn Gene IDs,
##    so a clone name whose record was merged into another gene resolves to
##    nothing above. The merge is taken from NCBI Gene and checked here only
##    as far as org.Hs.eg.db can: the surviving ID must still carry the
##    symbol we ship.
merged <- data.frame(
    original_symbol = "FLJ23867",
    withdrawn_id = "200058",
    current_id = "5768"
)
merged$shipped <- brs_genes$current_symbol[
    match(merged$original_symbol, brs_genes$original_symbol)]
merged$symbol_now <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db, merged$current_id, "SYMBOL", "ENTREZID"))
cat("\nClone names merged into another NCBI Gene record:\n")
cat(sprintf("  %-12s Gene ID %s -> %s (%s today)  shipped: %s\n",
            merged$original_symbol, merged$withdrawn_id, merged$current_id,
            merged$symbol_now, merged$shipped), sep = "")
stopifnot(identical(unname(merged$symbol_now), merged$shipped))

## 4. Shape invariants.
stopifnot(
    nrow(brs_genes) == 71L,
    !anyDuplicated(brs_genes$original_symbol),
    !anyNA(brs_genes$current_symbol),
    !anyDuplicated(brs_genes$current_symbol),
    identical(as.integer(table(brs_genes$block)), c(13L, 58L)),
    !is.unsorted(brs_genes$original_symbol[brs_genes$block == 1L]),
    !is.unsorted(brs_genes$original_symbol[brs_genes$block == 2L])
)
cat("\nShape invariants OK.\n")

record_session("check_symbols")
