#' @name plotCountsPerBiotype
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit acidplots::plotCountsPerBiotype
#'
#' @inheritParams plotCounts
#' @inheritParams acidplots::params
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



## Updated 2019-07-23.
`plotCountsPerBiotype,bcbioRNASeq` <-  # nolint
    function(object, normalized) {
        validObject(object)
        normalized <- match.arg(normalized)
        args <- .dynamicTrans(object = object, normalized = normalized)
        do.call(
            what = plotCountsPerBiotype,
            args = matchArgsToDoCall(
                args = args,
                removeFormals = "normalized"
            )
        )
    }

f1 <- formals(`plotCountsPerBiotype,bcbioRNASeq`)
f2 <- methodFormals(
    f = "plotCountsPerBiotype",
    signature = "SummarizedExperiment",
    package = "acidplots"
)
f2 <- f2[setdiff(names(f2), c(names(f1), "assay"))]
f <- c(f1, f2)
## Ensure TPM is set first.
f[["normalized"]] <- unique(c("tpm", normalizedCounts))
formals(`plotCountsPerBiotype,bcbioRNASeq`) <- f



#' @rdname plotCountsPerBiotype
#' @export
setMethod(
    f = "plotCountsPerBiotype",
    signature = signature("bcbioRNASeq"),
    definition = `plotCountsPerBiotype,bcbioRNASeq`
)



## Updated 2019-07-23.
`plotCountsPerBroadClass,bcbioRNASeq` <-  # nolint
    function(object, normalized) {
        validObject(object)
        normalized <- match.arg(normalized)
        args <- .dynamicTrans(
            object = object,
            normalized = normalized
        )
        do.call(
            what = plotCountsPerBroadClass,
            args = matchArgsToDoCall(
                args = args,
                removeFormals = "normalized"
            )
        )
    }

f1 <- formals(`plotCountsPerBroadClass,bcbioRNASeq`)
f2 <- methodFormals(
    f = "plotCountsPerBroadClass",
    signature = "SummarizedExperiment",
    package = "acidplots"
)
f2 <- f2[setdiff(names(f2), c(names(f1), "assay"))]
f <- c(f1, f2)
## Ensure TPM is set first.
f[["normalized"]] <- unique(c("tpm", normalizedCounts))
f[["trans"]] <- trans
formals(`plotCountsPerBroadClass,bcbioRNASeq`) <- f



#' @rdname plotCountsPerBiotype
#' @export
setMethod(
    f = "plotCountsPerBroadClass",
    signature = signature("bcbioRNASeq"),
    definition = `plotCountsPerBroadClass,bcbioRNASeq`
)
