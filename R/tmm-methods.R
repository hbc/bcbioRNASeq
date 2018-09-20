#' **T**rimmed **M**ean of **M**-Values (TMM) Normalization
#'
#' TMM normalization is recommended for RNA-seq data generally when the majority
#' of genes are not differentially expressed. We use this as a quality control
#' tool when plotting counts per gene.
#'
#' @note Only gene-level counts are supported.
#'
#' @name tmm
#' @family Data Functions
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#'
#' @return `matrix`.
#'
#' @examples
#' x <- tmm(bcb_small)
#' summary(x)
NULL



#' @rdname tmm
#' @export
setMethod(
    f = "tmm",
    signature = signature("matrix"),
    definition = function(object) {
        message("Applying trimmed mean of M-values (TMM) normalization")
        object %>%
            DGEList() %>%
            calcNormFactors() %>%
            cpm(normalized.lib.sizes = TRUE)
    }
)



#' @rdname tmm
#' @export
setMethod(
    f = "tmm",
    signature = signature("SummarizedExperiment"),
    definition = function(object, assay = 1L) {
        assert_is_scalar(assay)
        validObject(object)
        tmm(assays(object)[[assay]])
    }
)



#' @rdname tmm
#' @export
setMethod(
    f = "tmm",
    signature = signature("bcbioRNASeq"),
    definition = function(object) {
        validObject(object)
        .assertIsGeneLevel(object)
        tmm(assay(object))
    }
)
