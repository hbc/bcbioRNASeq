# TODO Unit test with rle -- need to handle non-finite values better.



#' @name plotGenderMarkers
#' @inherit basejump.plots::plotGenderMarkers
#' @author Michael Steinbaugh
#'
#' @inheritParams params
#' @inheritParams basejump.globals::params
#'
#' @return `ggplot`.
#'
#' @examples
#' data(bcb)
#' plotGenderMarkers(bcb)
NULL



#' @importFrom basejump.generics plotGenderMarkers
#' @aliases NULL
#' @export
basejump.generics::plotGenderMarkers



# bcbioRNASeq ==================================================================
plotGenderMarkers.bcbioRNASeq <-  # nolint
    function(object, normalized) {
        validObject(object)
        normalized <- match.arg(normalized)
        counts <- counts(object, normalized = normalized)
        # Ensure counts are log2 scale.
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
f1 <- formals(plotGenderMarkers.bcbioRNASeq)
f2 <- methodFormals(
    f = "plotGenderMarkers",
    signature = "SummarizedExperiment",
    package = "basejump.plots"
)
f2 <- f2[setdiff(
    x = names(f2),
    y = c(names(f1), "assay", "countsAxisLabel")
)]
f <- c(f1, f2)
f[["normalized"]] <- normalizedCounts
formals(plotGenderMarkers.bcbioRNASeq) <- f



#' @rdname plotGenderMarkers
#' @export
setMethod(
    f = "plotGenderMarkers",
    signature = signature("bcbioRNASeq"),
    definition = plotGenderMarkers.bcbioRNASeq
)
