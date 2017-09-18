#' Trimmed Mean of M-Values (TMM) Normalization
#'
#' TMM normalization is recommended for RNA-seq data generally when the majority
#' of genes are not differentially expressed. We use this as a quality control
#' tool when plotting counts per gene.
#'
#' @rdname tmm
#' @name tmm
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return [matrix].
#' @export
#'
#' @examples
#' data(bcb, dds)
#'
#' # bcbioRNADataSet
#' tmm(bcb) %>% head
#'
#' \dontrun{
#' # DESeqDataSet
#' tmm(dds) %>% head
#'
#' # matrix
#' assay(bcb) %>% tmm %>% head
#' }
NULL



# Constructors ====
.tmm <- function(object) {
    object %>%
        as.matrix %>%
        DGEList %>%
        calcNormFactors %>%
        cpm(normalized.lib.sizes = TRUE)
}



# Methods ====
#' @rdname tmm
#' @export
setMethod("tmm", "bcbioRNADataSet", function(object) {
    assays(object)[["tmm"]]
})



#' @rdname tmm
#' @export
setMethod("tmm", "DESeqDataSet", function(object) {
    assay(object) %>%
        .tmm
})



#' @rdname tmm
#' @export
setMethod("tmm", "matrix", .tmm)
