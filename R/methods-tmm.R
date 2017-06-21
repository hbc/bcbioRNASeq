#' Trimmed mean of M-values (TMM) normalization
#'
#' TMM normalization is recommended for RNA-seq data generally when the majority
#' of genes are not differentially expressed. We use this as a quality control
#' tool when plotting counts per gene.
#'
#' @rdname tmm
#' @docType methods
#'
#' @author Michael Steinbaugh
#'
#' @param object Object containing counts.
#'
#' @export
setMethod("tmm", "bcbioRnaDataSet", function(object) {
    assays(object)[["tmm"]]
})

#' @rdname tmm
#' @export
setMethod("tmm", "matrix", function(object) {
    .tmm(object)
})

#' @rdname tmm
#' @export
setMethod("tmm", "DESeqDataSet", function(object) {
    object %>%
        counts %>%
        .tmm
})
