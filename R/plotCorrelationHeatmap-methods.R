#' @name plotCorrelationHeatmap
#' @author Michael Steinbaugh
#' @importMethodsFrom acidplots plotCorrelationHeatmap
#' @inherit acidplots::plotCorrelationHeatmap
#' @note Updated 2019-09-15.
#'
#' @inheritParams plotCounts
#' @inheritParams acidroxygen::params
#' @param ... Passthrough to `SummarizedExperiment` method defined in acidplots.
#'   See [acidplots::plotCorrelationHeatmap()] for details.
#'
#' @examples
#' data(bcb)
#' plotCorrelationHeatmap(bcb, method = "pearson")
#' plotCorrelationHeatmap(bcb, method = "spearman")
NULL



#' @rdname plotCorrelationHeatmap
#' @name plotCorrelationHeatmap
#' @importFrom bioverbs plotCorrelationHeatmap
#' @export
NULL



## Updated 2019-09-15.
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

formals(`plotCorrelationHeatmap,bcbioRNASeq`)[["normalized"]] <- dt



#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    f = "plotCorrelationHeatmap",
    signature = signature("bcbioRNASeq"),
    definition = `plotCorrelationHeatmap,bcbioRNASeq`
)
