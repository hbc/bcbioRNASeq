#' @name plotCorrelationHeatmap
#' @author Michael Steinbaugh
#' @inherit AcidPlots::plotCorrelationHeatmap
#' @note Updated 2022-03-07.
#'
#' @inheritParams plotCounts
#' @inheritParams AcidRoxygen::params
#' @param ... Passthrough to `SummarizedExperiment` method defined in AcidPlots.
#' See `AcidPlots::plotCorrelationHeatmap()` for details.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' plotCorrelationHeatmap(bcb, method = "pearson")
#' plotCorrelationHeatmap(bcb, method = "spearman")
NULL



## Updated 2020-09-15.
`plotCorrelationHeatmap,bcbioRNASeq` <- # nolint
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

formals(`plotCorrelationHeatmap,bcbioRNASeq`)[["normalized"]] <- # nolint
    .dt



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    f = "plotCorrelationHeatmap",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotCorrelationHeatmap,bcbioRNASeq`
)
