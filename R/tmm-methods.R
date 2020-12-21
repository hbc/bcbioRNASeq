#' @name tmm
#' @author Michael Steinbaugh
#' @inherit AcidGenerics::tmm
#' @note Updated 2020-01-20.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' x <- tmm(bcb)
#' summary(x)
NULL



## Updated 2020-01-20.
`tmm,matrix` <-  # nolint
    function(object) {
        cli_alert("Applying trimmed mean of M-values (TMM) normalization.")
        object <- DGEList(object)
        object <- calcNormFactors(object, method = "TMM")
        object <- cpm(object, normalized.lib.sizes = TRUE)
        object
    }



#' @rdname tmm
#' @export
setMethod(
    f = "tmm",
    signature = signature("matrix"),
    definition = `tmm,matrix`
)



## Updated 2019-07-23.
`tmm,SummarizedExperiment` <-  # nolint
    function(object) {
        validObject(object)
        tmm(counts(object))
    }



#' @rdname tmm
#' @export
setMethod(
    f = "tmm",
    signature = signature("SummarizedExperiment"),
    definition = `tmm,SummarizedExperiment`
)
