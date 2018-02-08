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
#' @importFrom viridis scale_color_viridis
.plotDEGenePCA.DESeqResults <- function(
    object,
    counts,
    interestingGroups = "sampleName",
    lfc = 0L,
    color = viridis::scale_color_viridis(discrete = TRUE),
    label = FALSE,
    returnData = FALSE) {
    # Passthrough: counts, interestingGroups, color, label, returnData
    # lfc
    if (!(is.numeric(lfc) && length(lfc) == 1L)) {
        abort("`lfc` must be a numeric string")
    }

    # Get the DE gene vector using `resultsTables()`
    resTbl <- resultsTables(
        res,
        lfc = lfc,
        annotable = FALSE,
        summary = FALSE,
        write = FALSE)
    genes <- c(
        resTbl[["degLFCUp"]][["ensgene"]],
        resTbl[["degLFCDown"]][["ensgene"]]
    )

    .plotPCA.DESeqTransform(
        counts,
        interestingGroups = interestingGroups,
        genes = genes,
        color = color,
        label = label,
        returnData = returnData)
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
