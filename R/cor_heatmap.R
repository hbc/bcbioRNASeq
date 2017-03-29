#' Correlation matrix heatmap
#'
#' @author Michael Steinbaugh
#'
#' @import pheatmap
#' @import SummarizedExperiment
#' @importFrom stats cor
#'
#' @param bcbio bcbio list object
#' @param counts RNA-Seq counts
#' @param method Correlation coefficient (or covariance) to be computed
#' @param ... Optional passthrough to \code{pheatmap()}
#'
#' @return Heatmap
#' @export
cor_heatmap <- function(bcbio,
                        counts,
                        method = "pearson",
                        ...) {
    name <- deparse(substitute(counts))
    metadata <- import_metadata(bcbio)

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
                           annotation = metadata[, bcbio$intgroup],
                           show_colnames = FALSE,
                           ...)
}
