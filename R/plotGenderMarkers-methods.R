# FIXME Need to improve the formals to match basejump.



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
#' # bcbioRNASeq ====
#' object <- bcb_small
#' plotGenderMarkers(object, normalized = "vst")
#'
#' # DESeqTransform ====
#' object <- as(deseq_small, "DESeqTransform")
#' plotGenderMarkers(object)
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



.plotGenderMarkers.DESeqDataSet <-  # nolint
    function(object, ...) {
        validObject(object)
        counts <- log2(counts(object, normalized = TRUE) + 1L)
        countsAxisLabel <- "normalized counts (log2)"
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



.plotGenderMarkers.DESeqTransform <-  # nolint
    function(object, ...) {
        validObject(object)
        rse <- as(object, "RangedSummarizedExperiment")
        do.call(
            what = plotGenderMarkers,
            args = list(
                object = rse,
                countsAxisLabel = .transformCountsAxisLabel(object),
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



#' @rdname plotGenderMarkers
#' @export
setMethod(
    f = "plotGenderMarkers",
    signature = signature("DESeqDataSet"),
    definition = .plotGenderMarkers.DESeqDataSet
)



#' @rdname plotGenderMarkers
#' @export
setMethod(
    f = "plotGenderMarkers",
    signature = signature("DESeqTransform"),
    definition = .plotGenderMarkers.DESeqTransform
)
