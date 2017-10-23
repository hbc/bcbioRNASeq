#' Trimmed Mean of M-Values (TMM) Normalization
#'
#' TMM normalization is recommended for RNA-seq data generally when the majority
#' of genes are not differentially expressed. We use this as a quality control
#' tool when plotting counts per gene.
#'
#' @rdname tmm
#' @name tmm
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @inheritParams AllGenerics
#'
#' @return [matrix].
#' @export
#'
#' @examples
#' # bcbioRNASeq
#' tmm(bcb) %>% summary()
#'
#' # DESeqDataSet
#' \dontrun{
#' tmm(dds) %>% summary()
#' }
#'
#' # matrix
#' \dontrun{
#' assay(bcb) %>%
#'     tmm() %>%
#'     summary()
#' }
NULL



# Constructors ====
#' @importFrom edgeR calcNormFactors cpm DGEList
.tmm <- function(object) {
    object %>%
        as.matrix() %>%
        DGEList() %>%
        calcNormFactors() %>%
        cpm(normalized.lib.sizes = TRUE)
}



# Methods ====
#' @rdname tmm
#' @export
setMethod(
    "tmm",
    signature("bcbioRNASeq"),
    function(object) {
        assays(object)[["tmm"]]
    })



#' @rdname tmm
#' @export
setMethod(
    "tmm",
    signature("DESeqDataSet"),
    function(object) {
        assay(object) %>%
            .tmm()
    })



#' @rdname tmm
#' @export
setMethod(
    "tmm",
    signature("matrix"),
    .tmm)
