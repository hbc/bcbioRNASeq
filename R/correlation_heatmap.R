#' Correlation matrix heatmap
#'
#' @author Michael Steinbaugh
#'
#' @import DESeq2
#' @import pheatmap
#' @import SummarizedExperiment
#' @importFrom stats cor
#'
#' @param bcbio bcbio list object
#' @param counts RNA-Seq counts as \code{DESeqDataSet} or \code{DESeqTransform}
#'   object. If input as \code{DESeqDataSet}, an \code{rlog} transformation will
#'   be applied.
#' @param method Correlation coefficient (or covariance) to be computed.
#'   Defaults to \code{pearson} but \code{spearman} can also be used.
#' @param ... Optional passthrough to \code{pheatmap()}
#'
#' @return Heatmap
#' @export
correlation_heatmap <- function(bcbio,
                                counts,
                                method = "pearson",
                                ...) {
    check_bcbio_object(bcbio)
    if (class(counts)[1] == "DESeqDataSet") {
        counts <- DESeq2::rlog(counts)
    }
    if (class(counts)[1] == "DESeqTransform") {
        counts <- SummarizedExperiment::assay(counts)
    } else {
        stop("counts must be input as DESeqDataSet or DESeqTransform object")
    }
    if (!method %in% c("pearson", "spearman")) {
        stop("invalid correlation regression method.
             must use pearson or spearman.")
    }

    name <- deparse(substitute(counts))
    metadata <- import_metadata(bcbio)

    stats::cor(counts, method = method) %>%
        pheatmap::pheatmap(main = paste(method,
                                        "correlation:",
                                        name),
                           annotation = metadata[, bcbio$intgroup],
                           show_colnames = FALSE,
                           ...)
}
