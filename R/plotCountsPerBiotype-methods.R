#' @name plotCountsPerBiotype
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit bioverbs::plotCountsPerBiotype
#' @inheritParams basejump::params
#' @inheritParams params
#' @examples
#' data(bcb)
#' plotCountsPerBiotype(bcb)
NULL



#' @importFrom bioverbs plotCountsPerBiotype
#' @aliases NULL
#' @export
bioverbs::plotCountsPerBiotype

#' @importFrom bioverbs plotCountsPerBroadClass
#' @aliases NULL
#' @export
bioverbs::plotCountsPerBroadClass



plotCountsPerBiotype.bcbioRNASeq <-  # nolint
    function(
        object,
        normalized
    ) {
        validObject(object)
        normalized <- match.arg(normalized)

        # Coerce to RSE.
        rse <- as(object, "RangedSummarizedExperiment")
        counts <- counts(object, normalized = normalized)
        assays(rse) <- list(counts)
        assayNames(rse) <- normalized

        # Set the counts axis label.
        countsAxisLabel <- paste(normalized, "counts")
        trans <- .normalizedTrans(normalized)
        if (trans != "identity") {
            countsAxisLabel <- paste(trans, countsAxisLabel)
        }

        do.call(
            what = plotCountsPerBiotype,
            args = matchArgsToDoCall(
                args = list(
                    object = rse,
                    trans = trans,
                    countsAxisLabel = countsAxisLabel
                ),
                removeFormals = "normalized"
            )
        )
    }

f1 <- formals(plotCountsPerBiotype.bcbioRNASeq)
f2 <- methodFormals(
    f = "plotCountsPerBiotype",
    signature = "SummarizedExperiment",
    package = "basejump"
)
f2 <- f2[setdiff(names(f2), c(names(f1), "assay"))]
f <- c(f1, f2)
# Ensure TPM is set first.
f[["normalized"]] <- unique(c("tpm", normalizedCounts))
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
        normalized
    ) {
        validObject(object)
        normalized <- match.arg(normalized)

        # Coerce to RSE.
        rse <- as(object, "RangedSummarizedExperiment")
        counts <- counts(object, normalized = normalized)
        assays(rse) <- list(counts)
        assayNames(rse) <- normalized

        # Set the counts axis label.
        countsAxisLabel <- paste(normalized, "counts")
        trans <- .normalizedTrans(normalized)
        if (trans != "identity") {
            countsAxisLabel <- paste(trans, countsAxisLabel)
        }

        do.call(
            what = plotCountsPerBroadClass,
            args = matchArgsToDoCall(
                args = list(
                    object = rse,
                    trans = trans,
                    countsAxisLabel = countsAxisLabel
                ),
                removeFormals = "normalized"
            )
        )
    }

f1 <- formals(plotCountsPerBroadClass.bcbioRNASeq)
f2 <- methodFormals(
    f = "plotCountsPerBroadClass",
    signature = "SummarizedExperiment",
    package = "basejump"
)
f2 <- f2[setdiff(names(f2), c(names(f1), "assay"))]
f <- c(f1, f2)
# Ensure TPM is set first.
f[["normalized"]] <- unique(c("tpm", normalizedCounts))
formals(plotCountsPerBroadClass.bcbioRNASeq) <- f



#' @rdname plotCountsPerBiotype
#' @export
setMethod(
    f = "plotCountsPerBroadClass",
    signature = signature("bcbioRNASeq"),
    definition = plotCountsPerBroadClass.bcbioRNASeq
)
