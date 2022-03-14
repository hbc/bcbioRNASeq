#' @name relativeLogExpression
#' @author Lorena Pantano, Michael Steinbaugh
#' @inherit AcidGenerics::relativeLogExpression
#' @note Updated 2020-01-20.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' ## bcbioRNASeq ====
#' data(bcb)
#' relativeLogExpression(bcb)
NULL



## Updated 2020-01-20.
`relativeLogExpression,matrix` <- # nolint
    function(object) {
        alert("Applying relative log expression (RLE) normalization.")
        object <- DGEList(object)
        object <- calcNormFactors(object, method = "RLE")
        object <- cpm(object, normalized.lib.sizes = TRUE)
        object
    }



#' @rdname relativeLogExpression
#' @export
setMethod(
    f = "relativeLogExpression",
    signature = signature(object = "matrix"),
    definition = `relativeLogExpression,matrix`
)



## Updated 2019-07-23.
`relativeLogExpression,SummarizedExperiment` <- # nolint
    function(object) {
        relativeLogExpression(counts(object))
    }



#' @rdname relativeLogExpression
#' @export
setMethod(
    f = "relativeLogExpression",
    signature = signature(object = "SummarizedExperiment"),
    definition = `relativeLogExpression,SummarizedExperiment`
)
