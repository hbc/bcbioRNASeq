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
#' @param scale Character indicating if the values should be centered and scaled
#'   in either the `row` direction, `column` direction, or `none`.
#' @param annotationCol [data.frame] that specifies the annotations shown on the
#'   right side of the heatmap. Each row of this [data.frame] defines the
#'   features of the heatmap columns.
#' @param title Plot title.
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
#' # Genes as Ensembl identifiers
#' genes <- counts(bcb)[1L:50L, ] %>% rownames
#' plotGeneHeatmap(bcb, genes = genes)
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
    title = NULL,
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
            .[rownames(.) %in% genes, , drop = FALSE]
    }

    # Subset zero counts
    counts <- counts %>%
        .[rowSums(.) > 0L, , drop = FALSE]
    if (!is.matrix(counts) |
        nrow(counts) < 2L) {
        stop("Need at least 2 genes to cluster")
    }

    # Convert Ensembl gene identifiers to symbol names, if necessary
    if (nrow(counts) <= 100L) {
        showRownames <- TRUE
    } else {
        showRownames <- FALSE
    }
    if (isTRUE(showRownames)) {
        counts <- gene2symbol(counts)
    }

    # pheatmap will error if `NULL` title is passed as `main`
    if (is.null(title)) {
        title <- ""
    }

    pheatmap(counts,
             scale = scale,
             show_rownames = showRownames,
             annotation_col = annotationCol,
             main = title,
             ...)
}



# Methods ====
#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "bcbioRNADataSet", function(
    object, ...) {
    counts <- counts(object, normalized = "rlog")
    annotationCol <- colData(object) %>%
        # S4 DataFrame doesn't work with `pheatmap()`, so coerce to data.frame
        as.data.frame %>%
        .[, metadata(object)[["interestingGroups"]], drop = FALSE]
    .plotGeneHeatmap(counts, annotationCol = annotationCol, ...)
})



#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "DESeqDataSet", function(object, ...) {
    counts(object, normalized = TRUE) %>%
        .plotGeneHeatmap(...)
})



#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "DESeqTransform", function(object, ...) {
    assay(object) %>%
        .plotGeneHeatmap(...)
})



#' @rdname plotGeneHeatmap
#' @export
setMethod("plotGeneHeatmap", "matrix", .plotGeneHeatmap)
