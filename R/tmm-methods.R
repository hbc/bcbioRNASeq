#' @name tmm
#' @author Michael Steinbaugh
#' @inherit basejump::tmm
#' @inheritParams params
#' @examples
#' data(bcb)
#' x <- tmm(bcb)
#' summary(x)
NULL



#' @importFrom basejump tmm
#' @aliases NULL
#' @export
basejump::tmm



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
