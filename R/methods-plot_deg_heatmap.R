#' Differentially Expressed Gene Heatmap
#'
#' This function is a simplified version of [plot_gene_heatmap()] that is
#' optimized for handling a [DESeqResults] object rather a gene vector. All of
#' the optional parameters for [plot_gene_heatmap()] are also available to this
#' function.
#'
#' @rdname plot_deg_heatmap
#' @author Michael Steinbaugh
#' @family Heatmaps
#'
#' @param counts Secondary object containing normalized counts.
#' @param alpha Alpha level cutoff.
#' @param lfc [log2] fold change ratio cutoff.
#'
#' @return Graphical output only.
#'
#' @examples
#' data(bcb)
#' dds <- DESeqDataSetFromTximport(
#'     txi = txi(bcb),
#'     colData = colData(bcb),
#'     design = formula(~group)) %>%
#'     DESeq
#' res <- results(dds)
#' rld <- rlog(dds)
#' plot_deg_heatmap(res, rld)



# FIXME Draft support...needs update
#' @rdname plot_deg_heatmap
.plot_deg_heatmap <- function(
    object,
    counts,
    alpha = 0.05,
    lfc = 0L,
    ...) {
    alpha <- metadata(object)[["alpha"]]
    genes <- object %>%
        as.data.frame %>%
        rownames_to_column("ensgene") %>%
        snake %>%
        filter(.data[["padj"]] < !!alpha,
               .data[["log2_fold_change"]] > !!lfc |
                   .data[["log2_fold_change"]] < -UQ(lfc)) %>%
        pull("ensgene") %>%
        sort
    plot_gene_heatmap(counts, genes = genes, ...)
}



#' @rdname plot_deg_heatmap
#' @export
setMethod(
    "plot_deg_heatmap",
    signature(object = "DESeqResults",
              counts = "DESeqTransform"),
    function(
        object,
        counts) {
        print(class(object))
        print(class(counts))
    })
