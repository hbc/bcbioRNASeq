#' @name relativeLogExpression
#' @author Lorena Pantano, Michael Steinbaugh
#' @inherit bioverbs::relativeLogExpression
#'
#' @inheritParams params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' relativeLogExpression(bcb)
NULL



#' @rdname relativeLogExpression
#' @name relativeLogExpression
#' @importFrom bioverbs relativeLogExpression
#' @usage relativeLogExpression(object, ...)
#' @export
NULL



## Last modified 2019-06-06.
relativeLogExpression.matrix <-
    function(object) {
        message("Applying relative log expression (RLE) normalization.")
        object %>%
            DGEList() %>%
            calcNormFactors(method = "RLE") %>%
            cpm(normalized.lib.sizes = TRUE)
    }



#' @rdname relativeLogExpression
#' @export
setMethod(
    f = "relativeLogExpression",
    signature = signature("matrix"),
    definition = relativeLogExpression.matrix
)



relativeLogExpression.SummarizedExperiment <-  # nolint
    function(object) {
        relativeLogExpression(counts(object))
    }



#' @rdname relativeLogExpression
#' @export
setMethod(
    f = "relativeLogExpression",
    signature = signature("SummarizedExperiment"),
    definition = relativeLogExpression.SummarizedExperiment
)
