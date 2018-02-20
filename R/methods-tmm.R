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
#' @inheritParams general
#'
#' @return [matrix].
#' @export
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' tmm(bcb) %>% summary()
#'
#' # DESeqDataSet
#' dds <- bcbio(bcb, "DESeqDataSet")
#' tmm(dds) %>% summary()
#'
#' # matrix
#' counts <- counts(bcb)
#' tmm(counts) %>% summary()
NULL



# Constructors =================================================================
#' @importFrom edgeR calcNormFactors cpm DGEList
.tmm.matrix <- function(object) {  # nolint
    assert_is_matrix(object)
    inform("Performing trimmed mean of M-values (TMM) normalization")
    object %>%
        DGEList() %>%
        calcNormFactors() %>%
        cpm(normalized.lib.sizes = TRUE)
}



# Methods ======================================================================
#' @rdname tmm
#' @export
setMethod(
    "tmm",
    signature("bcbioRNASeq"),
    function(object) {
        data <- assays(object)[["tmm"]]
        assert_is_matrix(data)
        data
    })



#' @rdname tmm
#' @export
setMethod(
    "tmm",
    signature("DESeqDataSet"),
    function(object) {
        tmm(assay(object))
    })



#' @rdname tmm
#' @export
setMethod(
    "tmm",
    signature("matrix"),
    .tmm.matrix)
