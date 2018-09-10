#' **T**rimmed **M**ean of **M**-Values (TMM) Normalization
#'
#' TMM normalization is recommended for RNA-seq data generally when the majority
#' of genes are not differentially expressed. We use this as a quality control
#' tool when plotting counts per gene.
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
#' # bcbioRNASeq ====
#' tmm(bcb_small) %>% summary()
NULL



#' @rdname tmm
#' @export
setMethod(
    "tmm",
    signature("bcbioRNASeq"),
    function(object) {
        validObject(object)
        .assertIsGeneLevel(object)
        tmm(assay(object))
    }
)



#' @rdname tmm
#' @export
setMethod(
    "tmm",
    signature("DESeqDataSet"),
    function(object) {
        validObject(object)
        tmm(assay(object))
    }
)



#' @rdname tmm
#' @export
setMethod(
    "tmm",
    signature("matrix"),
    function(object) {
        message("Applying trimmed mean of M-values (TMM) normalization")
        object %>%
            DGEList() %>%
            calcNormFactors() %>%
            cpm(normalized.lib.sizes = TRUE)
    }
)
