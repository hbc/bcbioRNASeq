#' Trimmed mean of M-values (TMM) normalization
#'
#' TMM normalization is recommended for RNA-seq data generally when the majority
#' of genes are not differentially expressed. We use this as a quality control
#' tool when plotting counts per gene.
#'
#' @rdname tmm
#'
#' @author Michael Steinbaugh
#'
#' @param counts Counts matrix.
.tmm <- function(counts) {
    counts %>%
        DGEList %>%
        calcNormFactors %>%
        cpm(normalized.lib.sizes = TRUE)
}
