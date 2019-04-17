#' @name tmm
#' @author Michael Steinbaugh
#' @inherit bioverbs::tmm
#' @inheritParams params
#' @examples
#' data(bcb)
#' x <- tmm(bcb)
#' summary(x)
NULL



#' @rdname tmm
#' @name tmm
#' @importFrom bioverbs tmm
#' @export
NULL



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
