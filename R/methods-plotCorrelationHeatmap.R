#' Plot Correlation Heatmap
#'
#' This function calculates a correlation matrix based on feature expression per
#' sample.
#'
#' @name plotCorrelationHeatmap
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase plotCorrelationHeatmap
#'
#' @inherit bcbioBase::plotCorrelationHeatmap
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
#'
#' # Inferno palette
#' plotCorrelationHeatmap(
#'     bcb_small,
#'     color = inferno,
#'     legendColor = inferno
#' )
#'
#' # Default pheatmap palette
#' plotCorrelationHeatmap(
#'     bcb_small,
#'     color = NULL,
#'     legendColor = NULL
#' )
NULL



# Methods ======================================================================
#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    "plotCorrelationHeatmap",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = c("rlog", "vst", "tmm", "tpm"),
        ...
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        counts <- counts(object, normalized = normalized)

        # Coerce to RangedSummarizedExperiment
        rse <- as(object, "RangedSummarizedExperiment")
        assays(rse) <- list("counts" = counts)
        validObject(rse)

        plotCorrelationHeatmap(rse, ...)
    }
)



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    "plotCorrelationHeatmap",
    signature("DESeqDataSet"),
    function(
        object,
        normalized = TRUE,
        ...
    ) {
        validObject(object)
        counts <- counts(object, normalized = normalized)

        # Coerce to RangedSummarizedExperiment
        rse <- as(object, "RangedSummarizedExperiment")
        assays(rse) <- list("counts" = counts)
        validObject(rse)

        plotCorrelationHeatmap(rse, ...)
    }
)
