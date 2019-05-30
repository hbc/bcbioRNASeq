#' Sample PCA plot for transformed data
#'
#' Wrapper for [DESeq2::plotPCA()] that improves principal component analysis
#' (PCA) sample coloring and labeling.
#'
#' PCA (Jolliffe, et al., 2002) is a multivariate technique that allows us to
#' summarize the systematic patterns of variations in the data. PCA takes the
#' expression levels for genes and transforms it in principal component space,
#' reducing each sample into one point. Thereby, we can separate samples by
#' expression variation, and identify potential sample outliers. The PCA plot is
#' a way to look at how samples are clustering.
#'
#' @name plotPCA
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param ... Passthrough arguments to `SummarizedExperiment` method, defined
#'   in [acidplots::plotPCA()].
#'
#' @seealso
#' - [DESeq2::plotPCA()].
#' - `getMethod("plotPCA", "DESeqTransform")`
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioRNASeq ====
#' plotPCA(
#'     object = bcb_small,
#'     normalized = "vst",
#'     label = TRUE
#' )
#' plotPCA(
#'     object = bcb_small,
#'     normalized = "rlog",
#'     interestingGroups = "sampleName",
#'     label = FALSE
#' )
NULL



#' @rdname plotPCA
#' @name plotPCA
#' @importFrom BiocGenerics plotPCA
#' @usage plotPCA(object, ...)
#' @export
NULL



plotPCA.bcbioRNASeq <-  # nolint
    function(
        object,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
        ...
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        message(paste("Using", normalized, "counts"))
        rse <- as(object, "RangedSummarizedExperiment")
        assays(rse) <- list(counts(object, normalized = normalized))
        plotPCA(rse, ...)
    }



#' @rdname plotPCA
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("bcbioRNASeq"),
    definition = plotPCA.bcbioRNASeq
)
