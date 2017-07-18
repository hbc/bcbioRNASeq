#' Gene heatmap
#'
#' These functions facilitate heatmap plotting of a specified set of genes. By
#' default, row- and column-wise hierarchical clustering is performed using the
#' Ward method, but this behavior can be overrided by setting `cluster_rows` or
#' `cluster_cols` to `FALSE`. When column clustering is disabled, the columns
#' are sorted by the interesting groups (`interesting_groups`) specified in the
#' [bcbioRNADataSet] and then the sample names.
#'
#' @rdname plot_gene_heatmap
#' @family Heatmaps
#' @author Michael Steinbaugh
#'
#' @param object Object.
#' @param genes Character vector of specific gene identifiers to plot.
#' @param clustering_method Hierarchical clustering method. Accepts the same
#'   values as [stats::hclust()].
#' @param cluster_rows Use hierarchical clustering to arrange rows.
#' @param cluster_cols Use hierarchical clustering to arrange columns.
#' @param scale Character indicating if the values should be centered and scaled
#'   in either the `row` direction, `column` direction, or `none`.
#' @param interesting_groups *Optional*. Interesting groups to label with bars
#'   above heatmap. If `NULL`, defaults to `interesting_groups` defined in the
#'   [bcbioRNADataSet].
#' @param annotation *Optional*. Alternative annotation to use. Useful when
#'   plotting more than one column.
#' @param title *Optional*. Text to include in plot title.
#' @param ... Additional arguments.
#'
#' @seealso [pheatmap::pheatmap()].
#'
#' @return Graphical output only.
#' @export
setMethod("plot_gene_heatmap", "bcbioRNADataSet", function(
    object,
    genes,
    clustering_method = "ward.D2",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    interesting_groups = NULL,
    annotation = NULL,
    title = NULL,
    ...) {
    # Interesting groups
    if (is.null(interesting_groups)) {
        interesting_groups <- metadata(object)[["interesting_groups"]]
    }

    # Heatmap title (`main` parameter)
    if (!is.null(title)) {
        main <- title
    } else {
        main <- "gene heatmap"
    }

    # Transformed counts
    if (is.null(dt)) {
        # If NULL, use rlog counts from bcbioRNADataSet
        dt <- assays(object)[["rlog"]]
    }
    counts <- assay(dt) %>% .[genes, ]

    # `cluster_cols = FALSE`: Turn off sample hierarchical clustering and sort
    # by interesting groups then sample name, if preferred
    if (!isTRUE(cluster_cols)) {
        sorted_cols <- colData(object) %>%
            as.data.frame %>%
            arrange(!!!syms(c(interesting_groups, "description"))) %>%
            pull("description")
        counts <- counts[, sorted_cols]
    }

    # Change rownames to readable external gene names
    if (nrow(counts) <= 100L) {
        show_rownames <- TRUE
        gene2symbol <- gene2symbol(object)
        rownames(counts) <- gene2symbol[rownames(counts), "symbol"]
    } else {
        show_rownames <- FALSE
    }

    # Per sample annotations of interest
    if (is.null(annotation)) {
        annotation <- colData(object) %>%
            as.data.frame %>%
            .[colnames(counts), interesting_groups, drop = FALSE]
    }

    pheatmap(counts,
             annotation = annotation,
             cluster_cols = cluster_cols,
             cluster_rows = cluster_rows,
             clustering_method = clustering_method,
             main = main,
             scale = scale,
             show_colnames = TRUE,
             show_rownames = show_rownames,
             ...)
})
