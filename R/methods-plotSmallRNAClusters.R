#' Cluster Small RNA Samples
#'
#' @rdname plot_srna_clusters
#' @author Lorena Pantano, Michael Steinbaugh
#' @family Small RNA-Seq Utilities
#'
#' @param object Object.
#'
#' @return [ggplot].
#' @export
#'
#' @examples
#' data(bcb)
#' plot_srna_clusters(bcb)
setMethod("plot_srna_clusters", "bcbioRNADataSet", function(object) {
    counts <- assay(object) %>%
        # Require at least 3 counts
        .[rowSums(. > 0L) > 3L, ] %>%
        # DESeq requires integers
        round
    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = colData(object),
        design = ~1L)
    vst_counts <- rlog(dds) %>% assay
    annotation_col <- colData(object) %>%
        as.data.frame %>%
        .[, metadata(object)[["interesting_groups"]], drop = FALSE]
    pheatmap(vst_counts,
             annotation_col = annotation_col,
             clustering_distance_cols = "correlation",
             clustering_method = "ward.D",
             scale = "row",
             show_rownames = FALSE)
    plotPCA(object)
})
