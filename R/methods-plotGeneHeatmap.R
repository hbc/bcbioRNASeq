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
#'
#' # bcbioRNADataSet
#' plotGeneHeatmap(bcb)
#'
#' # Genes as symbols (default)
#' genes <- c("Rgs20", "Rab23", "Ncoa2", "Sulf1")
#' plotGeneHeatmap(bcb, genes = genes)
#'
#' # Genes as Ensembl identifiers
#' genes <- counts(bcb)[1L:50L, ] %>% rownames
#' plotGeneHeatmap(bcb, genes = genes, symbol = FALSE)
#'
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
    annotationCol = NULL,
    ...) {
    counts <- as.matrix(object)
    # Check for identifier mismatch. Do this before zero count subsetting.
    if (!is.null(genes)) {
        if (!all(genes %in% rownames(counts))) {
            stop(paste(
                "Genes missing from counts matrix:",
                toString(setdiff(genes, rownames(counts)))),
                call. = FALSE)
        }
        counts <- counts %>%
            .[rownames(.) %in% genes, ]
    }
    # Subset zero counts
    counts <- counts %>%
        .[rowSums(.) > 0L, ]
    if (!is.matrix(counts) |
        nrow(counts) < 2L) {
        stop("Need at least 2 genes to cluster")
    }
    if (nrow(counts) <= 100L) {
        showRownames <- TRUE
    } else {
        showRownames <- FALSE
    }
    pheatmap(counts,
             scale = scale,
             show_rownames = showRownames,
             annotation_col = annotationCol,
             ...)
}



# Methods ====
#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "bcbioRNADataSet", function(
    object, symbol = TRUE, ...) {
    counts <- counts(object, normalized = "rlog")
    if (isTRUE(symbol)) {
        counts <- gene2symbol(counts)
    }
    annotationCol <- colData(object) %>%
        # S4 DataFrame doesn't work with `pheatmap()`
        as.data.frame %>%
        .[, metadata(object)[["interestingGroups"]], drop = FALSE]
    .plotGeneHeatmap(counts, annotationCol = annotationCol, ...)
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
