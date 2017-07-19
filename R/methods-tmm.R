#' Trimmed Mean of M-Values (TMM) Normalization
#'
#' TMM normalization is recommended for RNA-seq data generally when the majority
#' of genes are not differentially expressed. We use this as a quality control
#' tool when plotting counts per gene.
#'
#' @rdname tmm
#' @author Michael Steinbaugh
#'
#' @return [matrix].
#' @export
#'
#' @examples
#' data(bcb)
#' tmm(bcb) %>% head
setMethod("tmm", "bcbioRNADataSet", function(object) {
    assays(object)[["tmm"]]
})



#' @rdname tmm
#' @usage NULL
.tmm <- function(object) {
    message("Generating TMM-normalized counts with edgeR")
    object %>%
        as.matrix %>%
        DGEList %>%
        calcNormFactors %>%
        cpm(normalized.lib.sizes = TRUE)
}



#' @rdname tmm
#' @export
setMethod("tmm", "DESeqDataSet", function(object) {
    object %>%
        counts %>%
        .tmm
})



#' @rdname tmm
#' @export
setMethod("tmm", "matrix", .tmm)
