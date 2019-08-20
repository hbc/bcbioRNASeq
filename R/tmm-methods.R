#' @name tmm
#' @author Michael Steinbaugh
#' @inherit bioverbs::tmm
#' @note Updated 2019-08-20.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' x <- tmm(bcb)
#' summary(x)
NULL



#' @rdname tmm
#' @name tmm
#' @importFrom bioverbs tmm
#' @usage tmm(object, ...)
#' @export
NULL



## Updated 2019-08-20.
`tmm,matrix` <-  # nolint
    function(object) {
        message("Applying trimmed mean of M-values (TMM) normalization.")
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
