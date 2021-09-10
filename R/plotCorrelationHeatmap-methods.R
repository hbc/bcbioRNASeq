#' @name plotCorrelationHeatmap
#' @author Michael Steinbaugh
#' @importMethodsFrom AcidPlots plotCorrelationHeatmap
#' @inherit AcidPlots::plotCorrelationHeatmap
#' @note Updated 2021-09-10.
#'
#' @inheritParams plotCounts
#' @inheritParams AcidRoxygen::params
#' @param ... Passthrough to `SummarizedExperiment` method defined in AcidPlots.
#'   See `AcidPlots::plotCorrelationHeatmap()` for details.
#'
#' @examples
#' ## bcbioRNASeq ====
#' data(bcb)
#' plotCorrelationHeatmap(bcb, method = "pearson")
#' plotCorrelationHeatmap(bcb, method = "spearman")
NULL



## Updated 2020-09-15.
`plotCorrelationHeatmap,bcbioRNASeq` <-  # nolint
    function(object, normalized, ...) {
        do.call(
            what = plotCorrelationHeatmap,
            args = list(
                object = .normalizedSE(
                    object = object,
                    normalized = match.arg(normalized)
                ),
                ...
            )
        )
    }

formals(`plotCorrelationHeatmap,bcbioRNASeq`)[["normalized"]] <- .dt



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    f = "plotCorrelationHeatmap",
    signature = signature("bcbioRNASeq"),
    definition = `plotCorrelationHeatmap,bcbioRNASeq`
)
