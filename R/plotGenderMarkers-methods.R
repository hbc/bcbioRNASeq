# FIXME
# Error in .local(object, ...) :
#     unused arguments (assay = 1, countsAxisLabel = "counts")
# Calls: plotGenderMarkers ... .local -> do.call -> do.call -> <Anonymous> -> <Anonymous>



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
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
        ...
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        counts <- counts(object, normalized = normalized)
        # Ensure counts are log2 scale
        if (!normalized %in% c("rlog", "vst")) {
            counts <- log2(counts + 1L)
        }
        countsAxisLabel <- paste(normalized, "counts (log2)")
        rse <- as(object, "RangedSummarizedExperiment")
        assay(rse) <- counts
        do.call(
            what = plotGenderMarkers,
            args = list(
                object = rse,
                countsAxisLabel = countsAxisLabel,
                ...
            )
        )
    }



#' @rdname plotGenderMarkers
#' @export
setMethod(
    f = "plotGenderMarkers",
    signature = signature("bcbioRNASeq"),
    definition = .plotGenderMarkers.bcbioRNASeq
)
