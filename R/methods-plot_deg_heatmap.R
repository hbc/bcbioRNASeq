#' Differentially expressed gene heatmap
#'
#' This function is a simplified version of [plot_gene_heatmap()] that is
#' optimized for handling a [DESeqResults] object rather a gene vector. All of
#' the optional parameters for [plot_gene_heatmap()] are also available to this
#' function.
#'
#' @rdname plot_deg_heatmap
#' @family Heatmaps
#' @author Michael Steinbaugh
#'
#' @param object Primary object containing a list of DEGs.
#' @param counts Secondary object containing normalized counts.
#' @param alpha Alpha level cutoff.
#' @param lfc [log2] fold change ratio cutoff.
#'
#' @export
setMethod("plot_deg_heatmap",
          signature(object = "bcbioRNAResults",
                    counts = "missingOrNULL"),
          function(
    object,
    counts,
    res,
    lfc,
    title = NULL,
    ...) {
    if (is.null(title)) {
        title <- "differentially expressed genes"
    }
    alpha <- metadata(res)[["alpha"]]
    genes <- res %>%
        as.data.frame %>%
        rownames_to_column("ensgene") %>%
        filter(.data[["padj"]] < alpha,
               .data[["log2FoldChange"]] > lfc |
                   .data[["log2FoldChange"]] < -lfc) %>%
        pull("ensgene") %>%
        sort
    plot_gene_heatmap(x, dt = dt, genes = genes, title = title, ...)
})



#' @rdname plot_deg_heatmap
#' @export
setMethod(
    "plot_deg_heatmap",
    signature(object = "DESeqResults", counts = "DESeqTransform"),
    function(
        object,
        counts) {
        # DESeqResults and DESeqTransform required
        print(class(object))
        print(class(counts))
    })
