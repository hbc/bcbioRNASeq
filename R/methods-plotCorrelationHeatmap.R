#' Plot Correlation Heatmap
#'
#' This function calculates a correlation matrix based on gene expression per
#' sample. By default, this function processes all gene counts per sample to
#' calculate the corrlation matrix. This behavior can be overrided with the
#' input of `gene` identifier vector. In this case, only the expression of the
#' desired genes will be used to calculate the correlation matrix.
#'
#' @name plotCorrelationHeatmap
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @importFrom basejump plotCorrelationHeatmap
#'
#' @inherit basejump::plotCorrelationHeatmap
#'
#' @inheritParams general
#'
#' @examples
#' # Pearson correlation
#' plotCorrelationHeatmap(bcb_small, method = "pearson")
#'
#' # Spearman correlation
#' plotCorrelationHeatmap(bcb_small, method = "spearman")
#'
#' # Inferno palette
#' plotCorrelationHeatmap(
#'     bcb_small,
#'     color = inferno,
#'     legendColor = inferno
#' )
#'
#' # Default pheatmap palette
#' plotCorrelationHeatmap(
#'     bcb_small,
#'     color = NULL,
#'     legendColor = NULL
#' )
NULL



# Methods ======================================================================
#' @rdname plotCorrelationHeatmap
#' @export
setMethod(
    "plotCorrelationHeatmap",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = c("rlog", "vst", "tmm", "tpm"),
        method = c("pearson", "spearman"),
        genes = NULL,
        samples = NULL,
        interestingGroups,
        color = viridis,
        legendColor = viridis,
        title = TRUE,
        ...
    ) {
        # Passthrough: method, genes, samples, color, legendColor
        normalized <- match.arg(normalized)
        method <- match.arg(method)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        assertFormalInterestingGroups(colData(object), interestingGroups)

        # Title
        if (isTRUE(title)) {
            title <- paste(
                normalized,
                paste(method, "correlation"),
                sep = " : "
            )
        }

        # Subset the counts matrix
        counts <- counts(object, normalized = normalized)
        if (is.character(genes)) {
            counts <- counts[genes, , drop = FALSE]
        }
        if (is.character(samples)) {
            counts <- counts[, samples, drop = FALSE]
            if (is.data.frame(annotationCol)) {
                annotationCol <- annotationCol[samples, , drop = FALSE]
            }
        }

        # Rename the matrix columns to `sampleName`, if necessary
        colData <- colData(object)
        if (!identical(
            x = colnames(object),
            y = as.character(colData[["sampleName"]])
        )) {
            rownames(colData) <- colData[["sampleName"]]
            colnames(counts) <- rownames(colData)
        }

        # Don't set annotation columns if we're only grouping by sample name
        if (identical(interestingGroups, "sampleName")) {
            annotationCol <- NULL
        } else {
            annotationCol <- colData(object) %>%
                .[, interestingGroups, drop = FALSE] %>%
                as.data.frame()
        }

        plotCorrelationHeatmap(
            object = counts,
            method = method,
            annotationCol = annotationCol,
            color = color,
            legendColor = legendColor,
            title = title,
            ...
        )
    }
)
