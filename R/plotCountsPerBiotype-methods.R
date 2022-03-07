#' @name plotCountsPerBiotype
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit AcidPlots::plotCountsPerBiotype
#' @note Updated 2022-03-07.
#'
#' @inheritParams plotCounts
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' plotCountsPerBiotype(bcb)
NULL



## Updated 2019-09-16.
`plotCountsPerBiotype,bcbioRNASeq` <-  # nolint
    function(object, normalized, ...) {
        do.call(
            what = plotCountsPerBiotype,
            args = .normalizedPlotArgs(
                object = object,
                normalized = match.arg(normalized),
                ...
            )
        )
    }

formals(`plotCountsPerBiotype,bcbioRNASeq`)[["normalized"]] <- .normalized



#' @rdname plotCountsPerBiotype
#' @export
setMethod(
    f = "plotCountsPerBiotype",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotCountsPerBiotype,bcbioRNASeq`
)



## Updated 2019-09-16.
`plotCountsPerBroadClass,bcbioRNASeq` <-  # nolint
    function(object, normalized, ...) {
        do.call(
            what = plotCountsPerBroadClass,
            args = .normalizedPlotArgs(
                object = object,
                normalized = match.arg(normalized),
                ...
            )
        )
    }

formals(`plotCountsPerBroadClass,bcbioRNASeq`) <-
    formals(`plotCountsPerBiotype,bcbioRNASeq`)



#' @rdname plotCountsPerBiotype
#' @export
setMethod(
    f = "plotCountsPerBroadClass",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotCountsPerBroadClass,bcbioRNASeq`
)
