#' @name plotCorrelationHeatmap
#' @inherit bioverbs::plotCorrelationHeatmap
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param ... Passthrough arguments to `SummarizedExperiment` method.
#'
#' @seealso
#' - `help("plotCorrelationHeatmap", "bcbioBase")`.
#' - `findMethod("plotCorrelationHeatmap", "SummarizedExperiment")`.
#'
#' @examples
#' # bcbioRNASeq ====
#' # Pearson correlation
#' plotCorrelationHeatmap(bcb_small, method = "pearson")
#'
#' # Spearman correlation
#' plotCorrelationHeatmap(bcb_small, method = "spearman")
NULL



#' @rdname plotCorrelationHeatmap
#' @name plotCorrelationHeatmap
#' @importFrom bioverbs plotCorrelationHeatmap
#' @usage plotCorrelationHeatmap(object, ...)
#' @export
NULL



plotCorrelationHeatmap.bcbioRNASeq <-  # nolint
    function(
        object,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
        ...
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts"))
        counts <- counts(object, normalized = normalized)

        # Coerce to RangedSummarizedExperiment
        rse <- as(object, "RangedSummarizedExperiment")
        assays(rse) <- list(counts = counts)
        validObject(rse)

        plotCorrelationHeatmap(rse, ...)
    }



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    f = "plotCorrelationHeatmap",
    signature = signature("bcbioRNASeq"),
    definition = plotCorrelationHeatmap.bcbioRNASeq
)



plotCorrelationHeatmap.DESeqTransform <-  # nolint
    function(object, ...) {
        validObject(object)
        counts <- assay(object)

        # Coerce to RangedSummarizedExperiment
        rse <- as(object, "RangedSummarizedExperiment")
        assays(rse) <- list(counts = counts)
        validObject(rse)

        plotCorrelationHeatmap(rse, ...)
    }



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    f = "plotCorrelationHeatmap",
    signature = signature("DESeqTransform"),
    definition = plotCorrelationHeatmap.DESeqTransform
)
