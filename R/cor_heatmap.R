#' Correlation matrix heatmap
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#' @import pheatmap
#' @import RColorBrewer
#' @import SummarizedExperiment
#' @importFrom grDevices colorRampPalette
#' @importFrom stats cor
#'
#' @param counts RNA-Seq counts
#' @param metadata Metadata data frame
#' @param factor Factor groups to use for clustering color headers
#' @param method Correlation coefficient (or covariance) to be computed
#' @param ... Passthrough to \code{pheatmap()}
#'
#' @return Heatmap
#' @export
cor_heatmap <- function(counts,
                        metadata,
                        factor = "group",
                        method = "pearson",
                        ...) {
    name <- deparse(substitute(counts))

    metadata <- as.data.frame(metadata)
    if (!is.data.frame(metadata)) {
        stop("A metadata data frame is required.")
    }
    if (!is.character(factor)) {
        stop("A factor group character vector is required.")
    }
    if (class(counts)[1] == "DESeqTransform") {
        counts <- SummarizedExperiment::assay(counts)
    }
    if (!is.matrix(counts)) {
        stop("A counts matrix is required.")
    }

    counts %>%
        stats::cor(., method = method) %>%
        pheatmap::pheatmap(main = paste(method,
                                        "correlation:",
                                        name),
                           annotation = metadata[, factor],
                           show_colnames = FALSE,
                           ...)
}
