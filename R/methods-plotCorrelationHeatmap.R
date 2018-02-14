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
#' @inherit plotHeatmap
#' @inheritParams counts
#' @inheritParams plotTotalReads
#'
#' @param method Correlation coefficient (or covariance) method to be computed.
#'   Defaults to `pearson` but `spearman` can also be used. Consult the
#'   [stats::cor()] documentation for more information.
#' @param samples *Optional.* Character vector of specific samples.
#' @param genes *Optional.* Character vector of specific gene identifiers to
#'   plot.
#'
#' @seealso
#' - [stats::cor()].
#' - [stats::hclust()].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # Pearson correlation
#' plotCorrelationHeatmap(bcb, method = "pearson")
#'
#' # Spearman correlation
#' plotCorrelationHeatmap(bcb, method = "spearman")
#'
#' # Inferno palette
#' plotCorrelationHeatmap(
#'     bcb,
#'     color = viridis::inferno(256),
#'     legendColor = viridis::inferno)
#'
#' # Default pheatmap palette
#' plotCorrelationHeatmap(
#'     bcb,
#'     color = NULL,
#'     legendColor = NULL)
NULL



# Constructors =================================================================
#' @importFrom dplyr mutate_all
#' @importFrom grDevices colorRampPalette
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats setNames
#' @importFrom S4Vectors cor
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom viridis viridis
.plotCorrelationHeatmap <- function(
    counts,
    method,
    annotationCol = NULL,
    genes = NULL,
    samples = NULL,
    color = viridis::viridis(256L),
    legendColor = viridis::viridis,
    title = TRUE) {
    # Check for supported correlation method
    validMethod <- c("pearson", "spearman")
    if (!method %in% validMethod) {
        abort(paste(
            "Supported methods:", toString(validMethod)
        ))
    }

    # Subset counts matrix by input genes, if desired
    if (!is.null(genes)) {
        counts <- counts[genes, , drop = FALSE]
    }

    # Subset count matrix by input samples, if desired
    if (!is.null(samples)) {
        counts <- counts[, samples, drop = FALSE]
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
    if (is.data.frame(annotationCol) & is.function(legendColor)) {
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

    # Set heatmap title (`main` parameter in pheatmap)
    if (isTRUE(title)) {
        title <- paste(method, "correlation")
    } else if (!is.character(title)) {
        title <- NULL
    }

    # If `color = NULL`, use the pheatmap default
    if (!is.character(color)) {
        color <- colorRampPalette(rev(
            brewer.pal(n = 7L, name = "RdYlBu")
        ))(100L)
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
            main = title,
            show_colnames = TRUE,
            show_rownames = TRUE)
}



# Methods ======================================================================
#' @rdname plotCorrelationHeatmap
#' @importFrom bcbioBase checkInterestingGroups
#' @importFrom viridis viridis
#' @export
setMethod(
    "plotCorrelationHeatmap",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = "rlog",
        method = "pearson",
        interestingGroups,
        genes = NULL,
        samples = NULL,
        color = viridis::viridis(256L),
        legendColor = viridis::viridis,
        title = TRUE) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        interestingGroups <- checkInterestingGroups(
            object = sampleMetadata(object),
            interestingGroups)

        counts <- counts(object, normalized = normalized)
        if (is.null(counts)) return(NULL)

        # Don't set annotation columns if we're only grouping by sample name
        if (identical(interestingGroups, "sampleName")) {
            annotationCol <- NULL
        } else {
            annotationCol <- colData(object) %>%
                .[, interestingGroups, drop = FALSE]
        }

        if (isTRUE(title) & is.character(normalized)) {
            title <- paste(normalized, method, "correlation")
        }

        .plotCorrelationHeatmap(
            counts = counts,
            method = method,
            annotationCol = annotationCol,
            genes = genes,
            samples = samples,
            color = color,
            legendColor = legendColor,
            title = title)
    })
