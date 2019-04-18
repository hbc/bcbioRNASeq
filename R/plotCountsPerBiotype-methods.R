#' @name plotCountsPerBiotype
#' @author Michael Steinbaugh, Rory Kirchner
#' @importMethodsFrom minimalism plotCountsPerBiotype
#' @inherit minimalism::plotCountsPerBiotype
#' @inheritParams minimalism::params
#' @inheritParams basejump::params
#' @inheritParams params
#' @examples
#' data(bcb)
#' plotCountsPerBiotype(bcb)
NULL



#' @rdname plotCountsPerBiotype
#' @name plotCountsPerBiotype
#' @importFrom bioverbs plotCountsPerBiotype
#' @export
NULL

#' @rdname plotCountsPerBiotype
#' @name plotCountsPerBroadClass
#' @importFrom bioverbs plotCountsPerBroadClass
#' @export
NULL



plotCountsPerBiotype.bcbioRNASeq <-  # nolint
    function(
        object,
        normalized,
        trans
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        trans <- match.arg(trans)
        args <- .dynamicTrans(
            object = object,
            normalized = normalized,
            trans = trans
        )
        do.call(
            what = plotCountsPerBiotype,
            args = matchArgsToDoCall(
                args = args,
                removeFormals = "normalized"
            )
        )
    }

f1 <- formals(plotCountsPerBiotype.bcbioRNASeq)
f2 <- methodFormals(
    f = "plotCountsPerBiotype",
    signature = "SummarizedExperiment",
    package = "minimalism"
)
f2 <- f2[setdiff(names(f2), c(names(f1), "assay"))]
f <- c(f1, f2)
# Ensure TPM is set first.
f[["normalized"]] <- unique(c("tpm", normalizedCounts))
f[["trans"]] <- trans
formals(plotCountsPerBiotype.bcbioRNASeq) <- f



#' @rdname plotCountsPerBiotype
#' @export
setMethod(
    f = "plotCountsPerBiotype",
    signature = signature("bcbioRNASeq"),
    definition = plotCountsPerBiotype.bcbioRNASeq
)



plotCountsPerBroadClass.bcbioRNASeq <-  # nolint
    function(
        object,
        normalized,
        trans
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        trans <- match.arg(trans)
        args <- .dynamicTrans(
            object = object,
            normalized = normalized,
            trans = trans
        )
        do.call(
            what = plotCountsPerBroadClass,
            args = matchArgsToDoCall(
                args = args,
                removeFormals = "normalized"
            )
        )
    }

f1 <- formals(plotCountsPerBroadClass.bcbioRNASeq)
f2 <- methodFormals(
    f = "plotCountsPerBroadClass",
    signature = "SummarizedExperiment",
    package = "minimalism"
)
f2 <- f2[setdiff(names(f2), c(names(f1), "assay"))]
f <- c(f1, f2)
# Ensure TPM is set first.
f[["normalized"]] <- unique(c("tpm", normalizedCounts))
f[["trans"]] <- trans
formals(plotCountsPerBroadClass.bcbioRNASeq) <- f



#' @rdname plotCountsPerBiotype
#' @export
setMethod(
    f = "plotCountsPerBroadClass",
    signature = signature("bcbioRNASeq"),
    definition = plotCountsPerBroadClass.bcbioRNASeq
)
