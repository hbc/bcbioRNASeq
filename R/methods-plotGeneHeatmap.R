#' Gene Heatmap
#'
#' These functions facilitate heatmap plotting of a specified set of genes. By
#' default, row- and column-wise hierarchical clustering is performed using the
#' Ward method, but this behavior can be overrided by setting `cluster_rows` or
#' `cluster_cols` to `FALSE`. When column clustering is disabled, the columns
#' are sorted by the interesting groups (`interestingGroups`) specified in the
#' [bcbioRNADataSet] and then the sample names.
#'
#' @rdname plotGeneHeatmap
#' @name plotGeneHeatmap
#'
#' @param object Counts matrix.
#' @param genes Character vector of specific gene identifiers to plot.
#' @param symbol Match against Ensembl gene symbols.
#' @param scale Character indicating if the values should be centered and scaled
#'   in either the `row` direction, `column` direction, or `none`.
#' @param ... Additional arguments, passed to [pheatmap()].
#'
#' @seealso [pheatmap::pheatmap()].
#'
#' @return Graphical output only.
#'
#' @examples
#' data(bcb, dds, rld)
#' genes <- counts(bcb) %>% rownames %>% .[1L:50L]
#'
#' # bcbioRNADataSet
#' plotGeneHeatmap(bcb)
#' plotGeneHeatmap(bcb, genes, symbol = FALSE)
#'
#' # DESeqDataSet
#' plotGeneHeatmap(dds)
#'
#' # DESeqTransform
#' plotGeneHeatmap(rld)
NULL



# Constructors ====
.plotGeneHeatmap <- function(
    object,
    genes = NULL,
    scale = "row",
    ...) {
    counts <- object %>%
        as.matrix %>%
        # Subset zero counts
        .[rowSums(.) > 0L, ]
    if (!is.null(genes)) {
        counts <- counts %>%
            .[rownames(.) %in% genes, ]
    }
    if (!is.matrix(counts) |
        nrow(counts) < 2L) {
        stop("Need at least 2 genes to cluster")
    }
    if (length(counts) == 0L) {
        return(NULL)
    }
    if (nrow(counts) <= 100L) {
        showRownames <- TRUE
    } else {
        showRownames <- FALSE
    }
    pheatmap(counts,
             scale = scale,
             show_rownames = showRownames,
             ...)
}



# Methods ====
#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "bcbioRNADataSet", function(
    object, ..., symbol = TRUE) {
    counts <- counts(object, normalized = "rlog")
    if (isTRUE(symbol)) {
        counts <- gene2symbol(counts)
    }
    annotationCol <- colData(object)
        .[, c("rowname", metadata(object)[["interestingGroups"]])]
    .plotGeneHeatmap(counts, annotation_col = annotationCol, ...)
})



#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "DESeqDataSet", function(object, ...) {
    counts <- counts(object, normalized = TRUE)
    .plotGeneHeatmap(counts, ...)
})



#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "DESeqTransform", function(object, ...) {
    counts <- assay(object)
    .plotGeneHeatmap(counts, ...)
})



#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "matrix", .plotGeneHeatmap)
