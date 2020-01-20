#' @name relativeLogExpression
#' @author Lorena Pantano, Michael Steinbaugh
#' @inherit acidgenerics::relativeLogExpression
#' @note Updated 2020-01-20.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' relativeLogExpression(bcb)
NULL



#' @rdname relativeLogExpression
#' @name relativeLogExpression
#' @importFrom acidgenerics relativeLogExpression
#' @usage relativeLogExpression(object, ...)
#' @export
NULL



## Updated 2020-01-20.
`relativeLogExpression,matrix` <-  # nolint
    function(object) {
        cli_alert("Applying relative log expression (RLE) normalization.")
        object <- DGEList(object)
        object <- calcNormFactors(object, method = "RLE")
        object <- cpm(object, normalized.lib.sizes = TRUE)
        object
    }



#' @rdname relativeLogExpression
#' @export
setMethod(
    f = "relativeLogExpression",
    signature = signature("matrix"),
    definition = `relativeLogExpression,matrix`
)



## Updated 2019-07-23.
`relativeLogExpression,SummarizedExperiment` <-  # nolint
    function(object) {
        relativeLogExpression(counts(object))
    }



#' @rdname relativeLogExpression
#' @export
setMethod(
    f = "relativeLogExpression",
    signature = signature("SummarizedExperiment"),
    definition = `relativeLogExpression,SummarizedExperiment`
)
