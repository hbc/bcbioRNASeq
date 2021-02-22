#' @name plotCountsPerFeature
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit AcidPlots::plotCountsPerFeature
#' @note Updated 2019-09-16.
#'
#' @inheritParams plotCounts
#' @inheritParams AcidRoxygen::params
#' @param ... Passthrough to `SummarizedExperiment` method defined in AcidPlots.
#'   See [AcidPlots::plotCountsPerFeature()] for details.
#'
#' @examples
#' data(bcb)
#' plotCountsPerFeature(bcb)
NULL



## Updated 2019-09-16.
`plotCountsPerFeature,bcbioRNASeq` <-  # nolint
    function(object, normalized, ...) {
        do.call(
            what = plotCountsPerFeature,
            args = .normalizedPlotArgs(
                object = object,
                normalized = match.arg(normalized),
                ...
            )
        )
    }

formals(`plotCountsPerFeature,bcbioRNASeq`)[["normalized"]] <-
    unique(c("tmm", .normalized))



#' @rdname plotCountsPerFeature
#' @export
setMethod(
    f = "plotCountsPerFeature",
    signature = signature("bcbioRNASeq"),
    definition = `plotCountsPerFeature,bcbioRNASeq`
)



## Note that this function is defined in F1000 v2 workflow paper.

#' @rdname plotCountsPerFeature
#' @export
plotCountDensity <- function(object, ...) {
    plotCountsPerFeature(object = object, geom = "density", ...)
}
