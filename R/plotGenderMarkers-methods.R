#' @name plotGenderMarkers
#' @inherit bioverbs::plotGenderMarkers
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @note Currently only *Homo sapiens* and *Mus musculus* genomes are supported.
#'
#' @inheritParams plotCounts
#' @inheritParams general
#'
#' @seealso [plotCounts()].
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioRNASeq ====
#' plotGenderMarkers(bcb_small, normalized = "vst")
#'
#' # DESeqTransform ====
#' vst_small <- DESeq2::varianceStabilizingTransformation(dds_small)
#' plotGenderMarkers(vst_small)
NULL



#' @rdname plotGenderMarkers
#' @name plotGenderMarkers
#' @importFrom bioverbs plotGenderMarkers
#' @usage plotGenderMarkers(object, ...)
#' @export
NULL



plotGenderMarkers.SummarizedExperiment <-  # nolint
    function(object, ...) {
        validObject(object)
        # Load the relevant internal gender markers data
        organism <- metadata(object)[["organism"]]
        assert_is_a_string(organism)
        markers <- bcbioRNASeq::gender_markers
        assert_is_subset(camel(organism), names(markers))
        markers <- markers[[camel(organism)]]

        # Get the X and Y chromosome marker genes
        xGenes <- markers %>%
            filter(!!sym("chromosome") == "X") %>%
            pull("geneID")
        yGenes <- markers %>%
            filter(!!sym("chromosome") == "Y") %>%
            pull("geneID")

        rse <- as(object, "RangedSummarizedExperiment")
        return <- "wide"
        xPlot <- plotCounts(
            object = rse,
            genes = xGenes,
            return = return,
            ...
        ) +
            ggtitle("X chromosome")
        yPlot <- plotCounts(
            object = rse,
            genes = yGenes,
            return = return,
            ...
        ) +
            ggtitle("Y chromosome")

        plotlist <- list(x = xPlot, y = yPlot)
        plot_grid(plotlist = plotlist)
    }



#' @rdname plotGenderMarkers
#' @export
setMethod(
    f = "plotGenderMarkers",
    signature = signature("SummarizedExperiment"),
    definition = plotGenderMarkers.SummarizedExperiment
)



plotGenderMarkers.bcbioRNASeq <-  # nolint
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
        plotGenderMarkers(
            object = rse,
            countsAxisLabel = countsAxisLabel,
            ...
        )
    }



#' @rdname plotGenderMarkers
#' @export
setMethod(
    f = "plotGenderMarkers",
    signature = signature("bcbioRNASeq"),
    definition = plotGenderMarkers.bcbioRNASeq
)



plotGenderMarkers.DESeqDataSet <-  # nolint
    function(object, ...) {
        validObject(object)
        counts <- log2(counts(object, normalized = TRUE) + 1L)
        countsAxisLabel <- "normalized counts (log2)"
        rse <- as(object, "RangedSummarizedExperiment")
        assay(rse) <- counts
        plotGenderMarkers(
            object = rse,
            countsAxisLabel = countsAxisLabel,
            ...
        )
    }



#' @rdname plotGenderMarkers
#' @export
setMethod(
    f = "plotGenderMarkers",
    signature = signature("DESeqDataSet"),
    definition = plotGenderMarkers.DESeqDataSet
)



plotGenderMarkers.DESeqTransform <-  # nolint
    function(object, ...) {
        validObject(object)
        if ("rlogIntercept" %in% colnames(mcols(object))) {
            normalized <- "rlog"
        } else {
            normalized <- "vst"
        }
        countsAxisLabel <- paste(normalized, "counts (log2)")
        rse <- as(object, "RangedSummarizedExperiment")
        plotGenderMarkers(
            object = rse,
            countsAxisLabel = countsAxisLabel,
            ...
        )
    }



#' @rdname plotGenderMarkers
#' @export
setMethod(
    f = "plotGenderMarkers",
    signature = signature("DESeqTransform"),
    definition = plotGenderMarkers.DESeqTransform
)
