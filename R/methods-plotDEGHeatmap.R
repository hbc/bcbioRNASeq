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
#' @param gene2symbol Convert Ensembl gene identifiers to symbol names.
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
    alpha = 0.01,
    lfc = 0L,
    ...) {
    res <- as.data.frame(object) %>%
        camel %>%
        # Keep genes that pass alpha cutoff
        .[!is.na(.[["padj"]]), , drop = FALSE] %>%
        .[.[["padj"]] < alpha, , drop = FALSE] %>%
        # Keep genes that pass log2 fold change cutoff
        .[!is.na(.[["log2FoldChange"]]), , drop = FALSE] %>%
        .[.[["log2FoldChange"]] > lfc |
              .[["log2FoldChange"]] < -lfc, , drop = FALSE]
    genes <- rownames(res)
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
