# FIXME Need to remove `assay` from formals.



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
#' plotGenderMarkers(bcb_small, normalized = "vst")
NULL



.plotGenderMarkers.bcbioRNASeq <-  # nolint
    function(
        object,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle")
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
f2 <- f2[setdiff(names(f2), names(f1))]
f <- c(f1, f2)
formals(.plotGenderMarkers.bcbioRNASeq) <- f



#' @rdname plotGenderMarkers
#' @export
setMethod(
    f = "plotGenderMarkers",
    signature = signature("bcbioRNASeq"),
    definition = .plotGenderMarkers.bcbioRNASeq
)
