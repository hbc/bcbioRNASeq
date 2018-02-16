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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
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
.tmm <- function(object) {
    assert_is_matrix(object)
    inform("Performing trimmed mean of M-values (TMM) normalization")
    object %>%
        DGEList() %>%
        calcNormFactors() %>%
        cpm(normalized.lib.sizes = TRUE)
}



.tmm.bcbioRNASeq <- function(object) {  # nolint
    tmm <- assays(object)[["tmm"]]
    assert_is_matrix(tmm, severity = "warning")
    if (!is.matrix(tmm)) {
        warn(paste(
            "TMM counts are not stashed in `assays()`.",
            "Recalculating on the fly.",
            "It is recomended that you restash these values in `assays()`."
        ))
        tmm <- tmm(counts(object, normalized = FALSE))
    }
    tmm
}



# Methods ======================================================================
#' @rdname tmm
#' @export
setMethod(
    "tmm",
    signature("bcbioRNASeq"),
    .tmm.bcbioRNASeq)



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
    .tmm)
