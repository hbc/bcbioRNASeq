#' Plot Correlation Heatmap
#'
#' This function calculates a correlation matrix based on gene expression per
#' sample. By default, this function processes all gene counts per sample to
#' calculate the corrlation matrix. This behavior can be overrided with the
#' input of `gene` identifier vector. In this case, only the expression of the
#' desired genes will be used to calculate the correlation matrix.
#'
#' @rdname plotCorrelationHeatmap
#' @author Michael Steinbaugh
#' @family Heatmaps
#'
#' @param transform String specifying `rlog` (**recommended**) or `vst`
#'   (`varianceStabilizingTransformation`) [DESeqTransform] object slotted
#'   inside the [bcbioRNADataSet].
#' @param method Correlation coefficient (or covariance) to be computed.
#'   Defaults to `pearson` but `spearman` can also be used. Consult the
#'   [stats::cor()] documentation for more information.
#' @param clusteringMethod Hierarchical clustering method. Accepts the same
#'   values as [stats::hclust()].
#' @param samples *Optional*. Character vector of specific samples.
#' @param genes *Optional*. Character vector of specific gene identifiers to
#'   plot.
#' @param interestingGroups *Optional*. Interesting groups to label with bars
#'   above heatmap. If `NULL`, defaults to `interestingGroups` defined in the
#'   [bcbioRNADataSet].
#' @param annotationCol *Optional*. Alternative annotation to use. Useful when
#'   plotting more than one column.
#' @param title *Optional*. Text to include in plot title.
#' @param ... Additional arguments, passed to [pheatmap::pheatmap()].
#'
#' @seealso
#' - [stats::cor()].
#' - [stats::hclust()].
#' - [pheatmap::pheatmap()].
#'
#' @return [pheatmap()].
#' @export
#'
#' @examples
#' data(bcb)
#' plotCorrelationHeatmap(bcb)
setMethod("plotCorrelationHeatmap", "bcbioRNADataSet", function(
    object,
    transform = "rlog",
    method = "pearson",
    clusteringMethod = "ward.D2",
    genes = NULL,
    samples = NULL,
    interestingGroups = NULL,
    annotationCol = NULL,
    title = NULL,
    ...) {
    # Check for supported correlation method
    if (!method %in% c("pearson", "spearman")) {
        stop("Supported methods: pearson, spearman")
    }

    # Interesting groups
    if (is.null(interestingGroups)) {
        interestingGroups <- metadata(object)[["interestingGroups"]]
    }

    # Per sample annotations of interest
    if (is.null(annotationCol)) {
        annotationCol <- colData(object) %>%
            .[, interestingGroups, drop = FALSE] %>%
            as.data.frame
    }

    # Set heatmap title (`main` parameter)
    if (!is.null(title)) {
        main <- title
    } else {
        main <- paste(method, "correlation")
    }

    # Transformed counts
    if (!transform %in% c("rlog", "vst")) {
        stop("DESeqTransform must be rlog or vst")
    }
    # Get count matrix from `assays` slot
    counts <- assays(object)[[transform]] %>% assay

    # Subset counts matrix by input genes, if desired
    if (!is.null(genes)) {
        counts <- counts[genes, ]
    }

    # Subset count matrix by input samples, if desired
    if (!is.null(samples)) {
        counts <- counts[, samples]
        annotationCol <- annotationCol[samples, ]
    }

    counts %>%
        cor(method = method) %>%
        pheatmap(
            annotation_col = annotationCol,
            clustering_method = clusteringMethod,
            clustering_distance_rows = "correlation",
            clustering_distance_cols = "correlation",
            main = main,
            show_colnames = FALSE,
            show_rownames = TRUE,
            ...)
})
