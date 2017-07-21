#' Gene Heatmap
#'
#' These functions facilitate heatmap plotting of a specified set of genes. By
#' default, row- and column-wise hierarchical clustering is performed using the
#' Ward method, but this behavior can be overrided by setting `cluster_rows` or
#' `cluster_cols` to `FALSE`. When column clustering is disabled, the columns
#' are sorted by the interesting groups (`interesting_groups`) specified in the
#' [bcbioRNADataSet] and then the sample names.
#'
#' @rdname plot_gene_heatmap
#' @author Michael Steinbaugh
#' @family Heatmaps
#'
#' @param counts Counts matrix.
#' @param genes Character vector of specific gene identifiers to plot.
#' @param symbol Match against Ensembl gene symbols.
#' @param cluster_rows Use hierarchical clustering to arrange rows.
#' @param cluster_cols Use hierarchical clustering to arrange columns.
#' @param scale Character indicating if the values should be centered and scaled
#'   in either the `row` direction, `column` direction, or `none`.
#' @param annotation *Optional*. Alternative annotation to use. Useful when
#'   plotting more than one column.
#' @param ... Additional arguments, passed to [pheatmap()].
#'
#' @seealso [pheatmap::pheatmap()].
#'
#' @return Graphical output only.
#'
#' @examples
#' data(bcb)
#' genes <- counts(bcb) %>% rownames %>% .[1L:50L]
#'
#' # bcbioRNADataSet
#' plot_gene_heatmap(bcb)
#' plot_gene_heatmap(bcb, genes, symbol = FALSE)
#'
#' # DESeqDataSet
#' dds <- DESeqDataSetFromTximport(
#'     txi = txi(bcb),
#'     colData = colData(bcb),
#'     design = formula(~group)) %>%
#'     DESeq
#' plot_gene_heatmap(dds)
#'
#' # DESeqTransform
#' rld <- rlog(dds)
#' plot_gene_heatmap(rld)



#' @rdname plot_gene_heatmap
.plot_gene_heatmap <- function(
    counts,
    genes = NULL,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    annotation = NULL,
    ...) {
    counts <- counts %>%
        as.matrix %>%
        # Subset zero counts
        .[rowSums(.) > 0, ]
    if (!is.null(genes)) {
        counts <- counts %>%
            .[rownames(.) %in% genes, ]
    }
    if (!is.matrix(counts) |
        nrow(counts) < 2L) {
        stop("Need at least 2 genes to cluster")
    }
    if (length(counts) == 0) return(NULL)
    if (nrow(counts) <= 100L) {
        show_rownames <- TRUE
    } else {
        show_rownames <- FALSE
    }
    pheatmap(counts,
             annotation = annotation,
             cluster_cols = cluster_cols,
             cluster_rows = cluster_rows,
             scale = scale,
             show_colnames = TRUE,
             show_rownames = show_rownames,
             ...)
}



#' @rdname plot_gene_heatmap
#' @export
setMethod("plot_gene_heatmap", "bcbioRNADataSet", function(
    object, ..., symbol = FALSE) {
    counts <- counts(object, normalized = "rlog")
    if (isTRUE(symbol)) {
        # Convert Ensembl gene identifiers to symbols
        rownames(counts) <- gene2symbol(object) %>%
            .[rownames(counts), "symbol"]
    }
    annotation <- colData(object) %>%
        as("tibble") %>%
        tidy_select(c("rowname", metadata(object)[["interesting_groups"]])) %>%
        as.data.frame %>%
        column_to_rownames
    .plot_gene_heatmap(counts, annotation = annotation, ...)
})



#' @rdname plot_gene_heatmap
#' @export
setMethod("plot_gene_heatmap", "DESeqDataSet", function(object, ...) {
    counts <- counts(object, normalized = TRUE)
    # FIXME Add annotation support based on colData
    .plot_gene_heatmap(counts, ...)
})



#' @rdname plot_gene_heatmap
#' @export
setMethod("plot_gene_heatmap", "DESeqTransform", function(object, ...) {
    counts <- assay(object)
    # FIXME Add annotation support based on colData
    .plot_gene_heatmap(counts, ...)
})
