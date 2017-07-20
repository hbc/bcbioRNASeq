#' Cluster Small RNA samples
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
    counts <- counts(object)
    # FIXME This slot is broken
    design <- metadata(object)[["metadata"]]
    if (is.null(design)) {
        warning("Broken code")
        return(NULL)
    }
    dds <- DESeqDataSetFromMatrix(
        counts[rowSums(counts > 0L) > 3L, ],
        colData = design,
        design = ~1L)
    vst <- rlog(dds, betaPriorVar = FALSE)
    annotation_col <- design %>%
        .[, metadata(object)[["interesting_groups"]], drop = FALSE]
    pheatmap(assay(vst),
             annotation_col = annotation_col,
             clustering_distance_cols = "correlation",
             clustering_method = "ward.D",
             scale = "row",
             show_rownames = FALSE)
    plot_pca(metadata(object), vst)
})
