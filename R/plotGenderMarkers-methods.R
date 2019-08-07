#' @name plotGenderMarkers
#' @author Michael Steinbaugh
#' @inherit acidplots::plotGenderMarkers
#'
#' @inheritParams plotCounts
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' plotGenderMarkers(bcb)
NULL



#' @rdname plotGenderMarkers
#' @name plotGenderMarkers
#' @importFrom bioverbs plotGenderMarkers
#' @importMethodsFrom acidplots plotGenderMarkers
#' @usage plotGenderMarkers(object, ...)
#' @export
NULL



## Updated 2019-07-23.
`plotGenderMarkers,bcbioRNASeq` <-  # nolint
    function(object, normalized) {
        validObject(object)
        normalized <- match.arg(normalized)
        counts <- counts(object, normalized = normalized)
        ## Ensure counts are log2 scale.
        if (!normalized %in% c("rlog", "vst")) {
            counts <- log2(counts + 1L)
        }
        countsAxisLabel <- paste(normalized, "counts (log2)")
        rse <- as(object, "RangedSummarizedExperiment")
        assay(rse) <- counts
        do.call(
            what = plotGenderMarkers,
            args = matchArgsToDoCall(
                args = list(
                    object = rse,
                    countsAxisLabel = countsAxisLabel
                ),
                removeFormals = "normalized"
            )
        )
    }

f1 <- formals(`plotGenderMarkers,bcbioRNASeq`)
f2 <- methodFormals(
    f = "plotGenderMarkers",
    signature = "SummarizedExperiment",
    package = "acidplots"
)
f2 <- f2[setdiff(
    x = names(f2),
    y = c(names(f1), "assay", "countsAxisLabel")
)]
f <- c(f1, f2)
f[["normalized"]] <- normalizedCounts
formals(`plotGenderMarkers,bcbioRNASeq`) <- f



#' @rdname plotGenderMarkers
#' @export
setMethod(
    f = "plotGenderMarkers",
    signature = signature("bcbioRNASeq"),
    definition = `plotGenderMarkers,bcbioRNASeq`
)
