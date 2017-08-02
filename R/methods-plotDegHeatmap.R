#' Differentially Expressed Gene Heatmap
#'
#' This function is a simplified version of [plotGeneHeatmap()] that is
#' optimized for handling a [DESeqResults] object rather a gene vector. All of
#' the optional parameters for [plotGeneHeatmap()] are also available to this
#' function.
#'
#' @rdname plotDEGHeatmap
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
#' plotDEGHeatmap(res, rld)



#' @rdname plotDEGHeatmap
.plotDEGHeatmap <- function(
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
               .data[["log2FoldChange"]] > !!lfc |
                   .data[["log2FoldChange"]] < -UQ(lfc)) %>%
        pull("ensgene") %>%
        sort
    plotGeneHeatmap(counts, genes = genes, ...)
}



#' @rdname plotDEGHeatmap
#' @export
setMethod(
    "plotDEGHeatmap",
    signature(object = "DESeqResults",
              counts = "DESeqTransform"),
    function(
        object,
        counts) {
        print(class(object))
        print(class(counts))
    })
