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
#' @param raw_counts Raw counts matrix
#'
#' @export
tmm_normalize <- function(raw_counts) {
    raw_counts %>%
        DGEList %>%
        calcNormFactors %>%
        cpm(normalized.lib.sizes = TRUE)
}
