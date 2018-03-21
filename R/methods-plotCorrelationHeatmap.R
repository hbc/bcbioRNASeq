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
#' @inherit plotHeatmap
#'
#' @inheritParams counts
#' @inheritParams plotTotalReads
#' @param method Correlation coefficient (or covariance) method to be computed.
#'   Defaults to "`pearson`" but "`spearman`" can also be used. Consult the
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
        interestingGroups,
        genes = NULL,
        samples = NULL,
        color = viridis,
        legendColor = viridis,
        title = TRUE,
        ...
    ) {
        # Passthrough: method, genes, samples, color, legendColor
        normalized <- match.arg(normalized)
        counts <- counts(object, normalized = normalized)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        assertFormalInterestingGroups(colData(object), interestingGroups)

        # Title
        if (isTRUE(title)) {
            title <- paste(normalized, method, "correlation")
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
            genes = genes,
            samples = samples,
            color = color,
            legendColor = legendColor,
            title = title,
            ...
        )
    }
)
