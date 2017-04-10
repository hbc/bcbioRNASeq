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
#'
#' @return Heatmap
#' @export
correlation_heatmap <- function(bcbio,
                                counts,
                                method = "pearson") {
    check_bcbio(bcbio)
    name <- deparse(substitute(counts))
    metadata <- import_metadata(bcbio)

    # Transform DESeqDataSet if necessary
    if (check_dds(counts, stop = FALSE)) {
        counts <- rlog(counts)
    }

    # Get counts from DESeqTransform object
    if (check_dt(counts)) {
        counts <- SummarizedExperiment::assay(counts)
    }

    # Pearson or Spearman correlation methods are supported
    if (!method %in% c("pearson", "spearman")) {
        stop("invalid correlation regression method.
             must use pearson or spearman.")
    }

    cor(counts, method = method) %>%
        pheatmap(main = paste(method,
                              "correlation:",
                              name),
                 annotation = metadata[, bcbio$intgroup],
                 show_colnames = TRUE,
                 show_rownames = TRUE)
}
