# Provenance: substantial parts of this file were written by Claude Opus 5
# (Anthropic) working from the published description of the method, under the
# direction of and reviewed by the package authors, who are responsible for
# its correctness and maintenance. See the Provenance section of README.md.
# Assisted-by: Claude Opus 5 and Claude Opus 5.5 (Anthropic)

#' BRAF-RAS Score (BRS) signature genes
#'
#' The 71-gene expression signature used to compute the BRAF-RAS Score (BRS)
#' for papillary thyroid carcinoma (PTC), as published in Agrawal et al.
#' (2014), Figure S7A. No machine-readable gene list was released with the
#' paper, so the symbols were transcribed from the published heatmap and then
#' verified against it row by row. The provenance and the checks that were
#' run are written up in `data-raw/brs_genes.md`, which is in the package
#' repository rather than the installed package:
#' <https://github.com/camilasauria/thyroidBRS/blob/main/data-raw/brs_genes.md>.
#'
#' @format A data frame with 71 rows and 4 columns:
#' \describe{
#'   \item{original_symbol}{Gene symbol exactly as it appears in the 2014
#'     publication's Figure S7A.}
#'   \item{current_symbol}{Current HGNC-approved symbol. Identical to
#'     `original_symbol` for genes whose symbol has not changed. Four were
#'     resolved against `org.Hs.eg.db`; `FLJ23867`, a cDNA clone name, was
#'     resolved from NCBI Gene to `QSOX1`. See the Alias resolution section.}
#'   \item{block}{Which of the two row-clusters of Figure S7A the gene belongs
#'     to: `1` for the 13-gene block, `2` for the 58-gene block. Each block is
#'     listed alphabetically in the figure, which is what makes the
#'     transcription checkable: no gene is out of order within its block.}
#'   \item{up_in}{Which reference group the gene is higher in, `"RAS"` for
#'     block 1 and `"BRAF"` for block 2. Refitting the centroids on TCGA-THCA
#'     recovers this split exactly: all 13 block-1 genes have a higher RAS
#'     centroid and all 58 block-2 genes a higher BRAF centroid, with
#'     no exceptions. That is an independent check on the gene list: a
#'     mis-transcribed symbol would land in the wrong block.}
#' }
#'
#' @section Alias resolution:
#' Five labels in the 2014 figure are not approved HGNC symbols today. Four
#' are former symbols: `FAM176A` (now `EVA1A`), `PVRL4` (`NECTIN4`), `TM7SF4`
#' (`DCSTAMP`) and `ARNTL` (`BMAL1`). Expression matrices built from current
#' GENCODE/HGNC annotation use the new symbols, so matching on
#' `original_symbol` silently loses those genes. Use `current_symbol`, which
#' is what [brs_fit()] does by default.
#'
#' `FLJ23867` is neither a symbol nor an alias in `org.Hs.eg.db`, HGNC or
#' GENCODE. It is a cDNA clone name: NCBI Gene record 200058
#' ("uncharacterized protein FLJ23867", 1q25.2) was withdrawn and merged into
#' 5768, `QSOX1` (quiescin sulfhydryl oxidase 1), whose record lists
#' FLJ23867 among its clone names. The merge dates from 2018, so in 2014
#' FLJ23867 was a Gene ID of its own and the published row may have measured
#' a transcript at that locus rather than the whole of `QSOX1`. `QSOX1` is
#' therefore the closest current gene, not necessarily the same feature. In
#' TCGA-THCA it is higher in BRAF-V600E than in RAS tumors and tracks the
#' rest of block 2.
#'
#' Resolution was checked against `org.Hs.eg.db`; re-run
#' `data-raw/check_symbols.R` when the annotation is updated. `data-raw/` is
#' not shipped in the installed package; it lives in the package repository,
#' <https://github.com/camilasauria/thyroidBRS>.
#'
#' @source Agrawal N, Akbani R, Aksoy BA, et al. "Integrated Genomic
#'   Characterization of Papillary Thyroid Carcinoma." Cell.
#'   2014;159(3):676-690. \doi{10.1016/j.cell.2014.09.050}, Figure S7A
#'   ("mRNA expression of the 71 genes used to derive the BRAF-V600E-RAS
#'   score (BRS) across 391 PTC samples").
#'
#' @examples
#' brs_genes
#' # Signature genes under current symbols:
#' brs_genes$current_symbol
#' # The two blocks of Figure S7A:
#' table(brs_genes$block)
#'
#' @export
brs_genes <- data.frame(
    original_symbol = c(
        # Block 1 of Figure S7A (13 genes, alphabetical)
        "ANKRD46", "CYB561", "GNA14", "HGD", "KATNAL2", "KCNAB1", "KCNIP3",
        "LGI3", "MLEC", "NQO1", "SFTPC", "SLC4A4", "SORBS2",
        # Block 2 of Figure S7A (58 genes, alphabetical)
        "ABTB2", "AHR", "ANKLE2", "ANXA1", "ANXA2P2", "ARNTL", "ASAP2",
        "BID", "CDC42EP1", "COL8A2", "CREB5", "CTSC", "CYP1B1", "DTX4",
        "DUSP5", "ETHE1", "FAM176A", "FAM20C", "FCHO1", "FLJ23867", "FN1",
        "FSTL3", "GABRB2", "GBP2", "ITGA3", "ITGB8", "KCNN4", "LAMB3",
        "LLGL1", "LY6E", "MDFIC", "MET", "PDE5A", "PDLIM4", "PLCD3",
        "PLEKHA6", "PNPLA5", "PPL", "PRICKLE1", "PROS1", "PTPRE", "PVRL4",
        "RASGEF1B", "RUNX1", "RUNX2", "SDC4", "SEL1L3", "SFTPB", "SLC35F2",
        "SOX4", "SPOCK2", "STK17B", "SYT12", "TACSTD2", "TBC1D2", "TGFBR1",
        "TM7SF4", "TMEM43"
    ),
    block = rep(c(1L, 2L), c(13L, 58L)),
    up_in = rep(c("RAS", "BRAF"), c(13L, 58L)),
    stringsAsFactors = FALSE
)

.brs_alias_to_current <- c(
    ARNTL = "BMAL1",
    FAM176A = "EVA1A",
    PVRL4 = "NECTIN4",
    TM7SF4 = "DCSTAMP",
    # Clone name, not an alias: NCBI Gene 200058 merged into 5768 (QSOX1).
    FLJ23867 = "QSOX1"
)

brs_genes$current_symbol <- ifelse(
    brs_genes$original_symbol %in% names(.brs_alias_to_current),
    .brs_alias_to_current[brs_genes$original_symbol],
    brs_genes$original_symbol
)

brs_genes <- brs_genes[, c(
    "original_symbol", "current_symbol", "block",
    "up_in"
)]
