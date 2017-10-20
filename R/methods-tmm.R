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
#' @importFrom edgeR calcNormFactors cpm DGEList
#'
#' @inheritParams AllGenerics
#'
#' @return [matrix].
#' @export
#'
#' @examples
#' data(bcb, dds)
#'
#' # bcbioRNASeq
#' tmm(bcb) %>%
#'     summary()
#'
#' \dontrun{
#' # DESeqDataSet
#' tmm(dds) %>%
#'     summary()
#'
#' # matrix
#' assay(bcb) %>%
#'     tmm() %>%
#'     summary()
#' }
NULL



# Constructors ====
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
    signature("bcbioRNASeqANY"),
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
