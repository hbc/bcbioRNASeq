# FIXME Need to remove `assay` from formals.
# FIXME Unit test with rle -- need to handle non-finite values better.



#' Plot Sexually Dimorphic Gender Marker Genes
#'
#' @name plotGenderMarkers
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#' @importFrom basejump plotGenderMarkers
#' @inherit basejump::plotGenderMarkers
#' @export
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' data(bcb_small)
#' plotGenderMarkers(bcb_small, normalized = "vst")
NULL



.plotGenderMarkers.bcbioRNASeq <-  # nolint
    function(
        object,
        normalized
    ) {
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
f1 <- formals(.plotGenderMarkers.bcbioRNASeq)
f2 <- methodFormals(
    f = "plotGenderMarkers",
    signature = "SummarizedExperiment"
)
f2 <- f2[setdiff(
    x = names(f2),
    y = c(names(f1), "assay", "countsAxisLabel")
)]
f <- c(f1, f2)
f[["normalized"]] <- normalizedCounts
formals(.plotGenderMarkers.bcbioRNASeq) <- f



#' @rdname plotGenderMarkers
#' @export
setMethod(
    f = "plotGenderMarkers",
    signature = signature("bcbioRNASeq"),
    definition = .plotGenderMarkers.bcbioRNASeq
)
