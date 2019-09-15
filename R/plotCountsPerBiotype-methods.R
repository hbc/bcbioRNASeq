#' @name plotCountsPerBiotype
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit acidplots::plotCountsPerBiotype
#' @note Updated 2019-09-15.
#'
#' @inheritParams plotCounts
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' plotCountsPerBiotype(bcb)
NULL



#' @rdname plotCountsPerBiotype
#' @name plotCountsPerBiotype
#' @importFrom bioverbs plotCountsPerBiotype
#' @importMethodsFrom acidplots plotCountsPerBiotype
#' @usage plotCountsPerBiotype(object, ...)
#' @export
NULL

#' @rdname plotCountsPerBiotype
#' @name plotCountsPerBroadClass
#' @importFrom bioverbs plotCountsPerBroadClass
#' @importMethodsFrom acidplots plotCountsPerBroadClass
#' @usage plotCountsPerBroadClass(object, ...)
#' @export
NULL



## Updated 2019-09-15.
`plotCountsPerBiotype,bcbioRNASeq` <-  # nolint
    function(object, normalized, ...) {
        do.call(
            what = plotCountsPerBiotype,
            args = .dynamicTrans(
                object = object,
                normalized = match.arg(normalized),
                ...
            )
        )
    }

f <- formals(`plotCountsPerBiotype,bcbioRNASeq`)
f[["normalized"]] <- unique(c("tpm", normalizedCounts))
formals(`plotCountsPerBiotype,bcbioRNASeq`) <- f



#' @rdname plotCountsPerBiotype
#' @export
setMethod(
    f = "plotCountsPerBiotype",
    signature = signature("bcbioRNASeq"),
    definition = `plotCountsPerBiotype,bcbioRNASeq`
)



## Updated 2019-09-15.
`plotCountsPerBroadClass,bcbioRNASeq` <-  # nolint
    function(object, normalized, ...) {
        do.call(
            what = plotCountsPerBroadClass,
            args = .dynamicTrans(
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
    signature = signature("bcbioRNASeq"),
    definition = `plotCountsPerBroadClass,bcbioRNASeq`
)
