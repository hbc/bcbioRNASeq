#' Differentially Expressed Gene Heatmap
#'
#' This function is a simplified version of [plotGeneHeatmap()] that is
#' optimized for handling a [DESeqResults] object rather a gene vector. All of
#' the optional parameters for [plotGeneHeatmap()] are also available to this
#' function.
#'
#' @rdname plotDEGHeatmap
#' @name plotDEGHeatmap
#'
#' @param counts Secondary object containing normalized counts.
#' @param alpha Alpha level cutoff.
#' @param lfc [log2] fold change ratio cutoff.
#'
#' @return Graphical output only.
#'
#' @examples
#' data(res, rld)
#' plotDEGHeatmap(res, rld)
NULL



# Constructors ====
.plotDEGHeatmap <- function(
    object,
    counts,
    alpha = 0.05,
    lfc = 0L,
    ...) {
    genes <- object %>%
        as.data.frame %>%
        rownames_to_column("ensgene") %>%
        camel %>%
        .[.[["padj"]] < alpha, ] %>%
        .[.[["log2FoldChange"]] > lfc |
              .[["log2FoldChange"]] < -lfc, ] %>%
        pull("ensgene") %>%
        sort
    plotGeneHeatmap(counts, genes = genes, ...)
}



# Methods ====
#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(object = "DESeqResults",
              counts = "DESeqTransform"),
    function(object, counts, ...) {
        alpha <- metadata(object)[["alpha"]]
        counts <- assay(counts)
        .plotDEGHeatmap(object, counts, alpha = alpha, ...)
    })



#' @rdname plotDEGHeatmap
#' @export
setMethod("plotDEGHeatmap", "data.frame", .plotDEGHeatmap)
