#' Trimmed Mean of M-Values (TMM) Normalization
#'
#' TMM normalization is recommended for RNA-seq data generally when the majority
#' of genes are not differentially expressed. We use this as a quality control
#' tool when plotting counts per gene.
#'
#' @name tmm
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `matrix`.
#' @export
#'
#' @examples
#' # bcbioRNASeq ====
#' tmm(bcb_small) %>% summary()
#'
#' # DESeqDataSet ====
#' tmm(dds_small) %>% summary()
#'
#' # matrix ====
#' counts(bcb_small) %>% tmm() %>% summary()
NULL



# Methods ======================================================================
#' @rdname tmm
#' @export
setMethod(
    "tmm",
    signature("bcbioRNASeq"),
    function(object) {
        tmm(assay(object))
    }
)



#' @rdname tmm
#' @export
setMethod(
    "tmm",
    signature("DESeqDataSet"),
    function(object) {
        tmm(assay(object))
    }
)



#' @rdname tmm
#' @importFrom edgeR calcNormFactors cpm DGEList
#' @export
setMethod(
    "tmm",
    signature("matrix"),
    function(object) {
        inform("Applying trimmed mean of M-values (TMM) normalization")
        object %>%
            DGEList() %>%
            calcNormFactors() %>%
            cpm(normalized.lib.sizes = TRUE)
    }
)
