#' Plot DEG PCA
#'
#' @rdname plotDEGenePCA
#' @name plotDEGenePCA
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @inheritParams plotPCA
#' @inheritParams plotDEGHeatmap
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "res.rda"),
#'     package = "bcbioRNASeq"))
#' load(system.file(
#'     file.path("extdata", "rld.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # DESeqResults + rlog counts
#' plotDEGenePCA(res, counts = rld)
NULL



# Constructors =================================================================
.plotDEGenePCA.DESeqResults <- function(
    object,
    counts,
    alpha,
    lfc = 0L,
    gene2symbol = TRUE,
    color = viridis::viridis(256L)) {
    if (missing(alpha)) {
        alpha <- metadata(object)[["alpha"]]
    }
    # FIXME resultsTables currently fails if annotable is FALSE
    resTbl <- resultsTables(
        res,
        annotable = FALSE,
        summary = FALSE,
        write = FALSE)
    genes <- c(
        resTbl[["degLFCUp"]][["ensgene"]],
        resTbl[["degLFCDown"]][["ensgene"]]
    )
    pcaCounts <- counts[genes, ]
    # DESeqTransform method
    plotPCA(counts, ntop = nrow(pcaCounts))
    # FIXME Use `returnData = TRUE`
    # stylize with the same settings as `plotPCA()`
}



# Methods ======================================================================
#' @rdname plotDEGenePCA
#' @export
setMethod(
    "plotDEGenePCA",
    signature(
        object = "DESeqResults",
        counts = "DESeqTransform"),
    .plotDEGenePCA.DESeqResults)
