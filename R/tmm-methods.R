#' Trimmed Mean of M-Values
#'
#' TMM normalization is recommended for RNA-seq data generally when the majority
#' of genes are not differentially expressed.
#'
#' @note Only recommended for gene-level counts.
#'
#' @name tmm
#' @author Michael Steinbaugh
#'
#' @inheritParams basejump.globals::params
#'
#' @return `matrix`.
#'
#' @references Robinson and Oshlack (2010).
#'
#' @seealso
#' - [edgeR::calcNormFactors()].
#' - [edgeR::cpm()].
#'
#' @examples
#' data(bcb_small)
#' x <- tmm(bcb_small)
#' summary(x)
NULL



# matrix =======================================================================
tmm.matrix <-  # nolint
    function(object) {
        message("Applying trimmed mean of M-values (TMM) normalization.")
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



# SummarizedExperiment =========================================================
tmm.SummarizedExperiment <-  # nolint
    function(object) {
        validObject(object)
        tmm(counts(object))
    }



#' @rdname tmm
#' @export
setMethod(
    f = "tmm",
    signature = signature("SummarizedExperiment"),
    definition = tmm.SummarizedExperiment
)
