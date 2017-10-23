#' Plot Correlation Heatmap
#'
#' This function calculates a correlation matrix based on gene expression per
#' sample. By default, this function processes all gene counts per sample to
#' calculate the corrlation matrix. This behavior can be overrided with the
#' input of `gene` identifier vector. In this case, only the expression of the
#' desired genes will be used to calculate the correlation matrix.
#'
#' @rdname plotCorrelationHeatmap
#' @name plotCorrelationHeatmap
#' @family Heatmaps
#' @author Michael Steinbaugh
#'
#' @inherit plotGeneHeatmap
#'
#' @param transform String specifying `rlog` (**recommended**) or `vst`
#'   (`varianceStabilizingTransformation`) [DESeqTransform] object slotted
#'   inside the [bcbioRNASeq] object.
#' @param method Correlation coefficient (or covariance) method to be computed.
#'   Defaults to `pearson` but `spearman` can also be used. Consult the
#'   [stats::cor()] documentation for more information.
#' @param samples *Optional*. Character vector of specific samples.
#' @param genes *Optional*. Character vector of specific gene identifiers to
#'   plot.
#'
#' @seealso
#' - [stats::cor()].
#' - [stats::hclust()].
#'
#' @examples
#' # Pearson correlation (default)
#' plotCorrelationHeatmap(bcb)
#'
#' # Spearman correlation
#' plotCorrelationHeatmap(bcb, method = "spearman")
#'
#' # Flip the palettes used for plot and legend
#' \dontrun{
#' plotCorrelationHeatmap(
#'     bcb,
#'     color = viridis::viridis(256),
#'     legendColor = viridis::inferno)
#' }
#'
#' # Default pheatmap palette
#' \dontrun{
#' plotCorrelationHeatmap(
#'     bcb,
#'     color = NULL,
#'     legendColor = NULL)
#' }
NULL



# Constructors ====
#' @importFrom dplyr mutate_all
#' @importFrom grDevices colorRampPalette
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats setNames
#' @importFrom S4Vectors cor
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom viridis inferno viridis
.plotCorrelationHeatmap <- function(
    counts,
    method,
    annotationCol = NULL,
    genes = NULL,
    samples = NULL,
    title = NULL,
    color = viridis::inferno(256),
    legendColor = viridis::viridis) {
    # Check for supported correlation method
    if (!method %in% c("pearson", "spearman")) {
        stop("Supported methods: pearson, spearman")
    }

    # Subset counts matrix by input genes, if desired
    if (!is.null(genes)) {
        counts <- counts[genes, , drop = FALSE]
    }

    # Subset count matrix by input samples, if desired
    if (!is.null(samples)) {
        counts <- counts[, samples]
        if (!is.null(annotationCol)) {
            annotationCol <- annotationCol[samples, , drop = FALSE]
        }
    }

    # Prepare the annotation columns
    if (!is.null(annotationCol)) {
        # Coerce annotation columns to factors
        annotationCol <- annotationCol %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            mutate_all(factor) %>%
            column_to_rownames()
    }

    # Define colors for each annotation column, if desired
    if (!is.null(annotationCol) & !is.null(legendColor)) {
        annotationColors <- lapply(
            seq_along(colnames(annotationCol)), function(a) {
                col <- annotationCol[[a]] %>%
                    levels()
                colors <- annotationCol[[a]] %>%
                    levels() %>%
                    length() %>%
                    legendColor
                names(colors) <- col
                colors
            }) %>%
            setNames(colnames(annotationCol))
    } else {
        annotationColors <- NULL
    }

    # Set heatmap title (`main` parameter)
    if (!is.null(title)) {
        main <- title
    } else {
        main <- paste(method, "correlation")
    }

    # If `color = NULL`, use the pheatmap default
    if (is.null(color)) {
        color <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    }

    counts %>%
        cor(method = method) %>%
        pheatmap(
            annotation_col = annotationCol,
            annotation_colors = annotationColors,
            border_color = NA,
            clustering_method = "ward.D2",
            clustering_distance_rows = "correlation",
            clustering_distance_cols = "correlation",
            color = color,
            main = main,
            show_colnames = FALSE,
            show_rownames = TRUE)
}



# Methods ====
#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    "plotCorrelationHeatmap",
    signature("bcbioRNASeq"),
    function(
        object,
        transform = "rlog",
        method = "pearson",
        genes = NULL,
        samples = NULL,
        title = NULL,
        color = viridis::inferno(256),
        legendColor = viridis::viridis) {
        # Transformed counts
        if (!transform %in% c("rlog", "vst")) {
            stop("DESeqTransform must be rlog or vst", call. = FALSE)
        }
        # Get count matrix from `assays` slot
        counts <- assays(object) %>%
            .[[transform]] %>%
            assay()
        interestingGroups <- interestingGroups(object)
        annotationCol <- colData(object) %>%
            .[, interestingGroups, drop = FALSE] %>%
            as.data.frame()
        .plotCorrelationHeatmap(
            counts = counts,
            method = method,
            annotationCol = annotationCol,
            genes = genes,
            samples = samples,
            title = title,
            color = color,
            legendColor = legendColor)
    })
