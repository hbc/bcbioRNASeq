#' Plot DEG PCA
#'
#' @rdname plotDEGenePCA
#' @name plotDEGenePCA
#' @author Michael Steinbaugh
#'
#' @inheritParams general
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



# Methods ======================================================================
#' @rdname plotDEGenePCA
#' @export
setMethod(
    "plotDEGenePCA",
    signature(
        object = "DESeqResults",
        counts = "DESeqTransform"),
    function(
        object,
        counts,
        interestingGroups = "sampleName",
        lfc = 0L,
        color = scale_color_viridis(discrete = TRUE),
        label = FALSE,
        returnData = FALSE) {
        # Passthrough: interestingGroups, color, label, returnData
        assert_is_a_number(lfc)
        assert_all_are_non_negative(lfc)

        # Get the DE gene vector using `resultsTables()`
        list <- resultsTables(
            object,
            lfc = lfc,
            summary = FALSE,
            write = FALSE)
        genes <- c(
            list[["degLFCUp"]][["ensgene"]],
            list[["degLFCDown"]][["ensgene"]])

        .plotPCA.DESeqTransform(
            object = counts,
            interestingGroups = interestingGroups,
            genes = genes,
            color = color,
            label = label,
            returnData = returnData)
    })



# Aliases ======================================================================
#' @rdname plotDEGenePCA
#' @export
plotDEGPCA <- function(...) {
    plotDEGenePCA(...)
}
