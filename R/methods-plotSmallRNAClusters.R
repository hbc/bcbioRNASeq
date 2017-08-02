#' Cluster Small RNA Samples
#'
#' @rdname plotSmallRNAClusters
#' @name plotSmallRNAClusters
#'
#' @return [ggplot].
#'
#' @examples
#' data(bcb)
#' plotSmallRNAClusters(bcb)
NULL



# Methods ====
#' @rdname plotSmallRNAClusters
#' @export
setMethod("plotSmallRNAClusters", "bcbioRNADataSet", function(object) {
    counts <- assay(object) %>%
        # Require at least 3 counts
        .[rowSums(. > 0L) > 3L, ] %>%
        # DESeq requires integers
        round
    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = colData(object),
        design = ~1L)
    vstCounts <- rlog(dds) %>% assay
    annotationCol <- colData(object) %>%
        as.data.frame %>%
        .[, metadata(object)[["interestingGroups"]], drop = FALSE]
    pheatmap(vstCounts,
             annotation_col = annotationCol,
             clustering_distance_cols = "correlation",
             clustering_method = "ward.D",
             scale = "row",
             show_rownames = FALSE)
    plotPCA(object)
})
