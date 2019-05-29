#' @name tmm
#' @inherit bioverbs::tmm
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param ... Additional arguments.
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



#' @rdname tmm
#' @name tmm
#' @importFrom bioverbs tmm
#' @usage tmm(object, ...)
#' @export
NULL



tmm.bcbioRNASeq <-  # nolint
    function(object) {
        validObject(object)
        tmm(assay(object))
    }



#' @rdname tmm
#' @export
setMethod(
    f = "tmm",
    signature = signature("bcbioRNASeq"),
    definition = tmm.bcbioRNASeq
)



tmm.DESeqDataSet <-  # nolint
    function(object) {
        validObject(object)
        tmm(assay(object))
    }



#' @rdname tmm
#' @export
setMethod(
    f = "tmm",
    signature = signature("DESeqDataSet"),
    definition = tmm.DESeqDataSet
)



tmm.matrix <-  # nolint
    function(object) {
        message("Applying trimmed mean of M-values (TMM) normalization")
        object %>%
            DGEList() %>%
            calcNormFactors() %>%
            cpm(normalized.lib.sizes = TRUE)
    }



#' @rdname tmm
#' @export
setMethod(
    f = "tmm",
    signature = signature("matrix"),
    definition = tmm.matrix
)
