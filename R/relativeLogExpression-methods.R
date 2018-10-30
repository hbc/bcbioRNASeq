# Don't attempt to abbreviate as `rle()`, since that is already defined in R.



#' Relative Log Expression
#'
#' @name relativeLogExpression
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @seealso [edgeR::calcNormFactors()].
#'
#' @references Anders and Huber (2010).
#'
#' @examples
#' data(bcb)
#' relativeLogExpression(bcb)
NULL



# matrix =======================================================================
relativeLogExpression.matrix <-  # nolint
    function(object) {
        t(t(object) / colMedians(object))
    }



#' @rdname relativeLogExpression
#' @export
setMethod(
    f = "relativeLogExpression",
    signature = signature("matrix"),
    definition = relativeLogExpression.matrix
)



# SummarizedExperiment =========================================================
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
