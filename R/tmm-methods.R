# TODO Move these methods to basejump.experiment.



#' @name tmm
#' @inherit basejump.generics::tmm
#' @author Michael Steinbaugh
#'
#' @inheritParams params
#'
#' @examples
#' data(bcb)
#' x <- tmm(bcb)
#' summary(x)
NULL



#' @importFrom basejump.generics tmm
#' @export
basejump.generics::tmm



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
