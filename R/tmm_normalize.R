#' Trimmed mean of M-values (TMM) normalization
#'
#' TMM normalization is recommended for RNA-Seq data generally when the majority
#' of genes are not differentially expressed. We use this as a quality control
#' tool when plotting counts per gene.
#'
#' @author Michael Steinbaugh
#'
#' @importFrom edgeR calcNormFactors cpm DGEList
#'
#' @param counts Counts matrix
#'
#' @export
tmm_normalize <- function(counts) {
    counts %>%
        edgeR::DGEList(.) %>%
        edgeR::calcNormFactors(.) %>%
        edgeR::cpm(normalized.lib.sizes = TRUE)
}
