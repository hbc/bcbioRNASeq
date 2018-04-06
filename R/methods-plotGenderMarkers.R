#' Plot Sexually Dimorphic Gender Markers
#'
#' @name plotGenderMarkers
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams plotGene
#' @inheritParams general
#' @param organism Organism name. Typically can be left unset and should be
#'   detected automatically, unless a spike-in FASTA sequence is provided
#'   containing a gene identifier that is first alphabetically in the count
#'   matrix rownames.
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioRNASeq ====
#' plotGenderMarkers(bcb_small)
#'
#' # DESeqDataSet ====
#' plotGenderMarkers(dds_small)
#'
#' # DESeqTransform ====
#' plotGenderMarkers(rld_small)
NULL



# Methods ======================================================================
#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("SummarizedExperiment"),
    function(object, ...) {
        validObject(object)
        return <- "wide"

        # Load the relevant internal gender markers data
        organism <- metadata(object)[["organism"]]
        if (!is_a_string(organism)) {
            organism <- detectOrganism(counts)
        }
        markers <- bcbioRNASeq::genderMarkers
        assert_is_subset(camel(organism), names(markers))
        markers <- markers[[camel(organism)]]

        xGenes <- markers %>%
            .[.[["chromosome"]] == "X", "geneID", drop = TRUE]
        yGenes <- markers %>%
            .[.[["chromosome"]] == "Y", "geneID", drop = TRUE]

        xPlot <- plotGene(
            object = rse,
            genes = xGenes,
            return = return,
            ...
        ) +
            ggtitle("X chromosome")

        yPlot <- plotGene(
            object = rse,
            genes = yGenes,
            return = return,
            ...
        ) +
            ggtitle("Y chromosome")

        plotlist <- list(xPlot, yPlot)
        plot_grid(plotlist = plotlist)
    }
)



#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = c("rlog", "vst", "tpm"),
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
)



#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("DESeqDataSet"),
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
)



#' @rdname plotGenderMarkers
#' @export
setMethod(
    "plotGenderMarkers",
    signature("DESeqTransform"),
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
)
