#' Relative Log Expression
#'
#' @name relativeLogExpression
#' @author Lorena Pantano, Michael Steinbaugh
#' @export
#'
#' @seealso [edgeR::calcNormFactors].
#'
#' @references Anders and Huber (2010).
#'
#' @examples
#' data(bcb_small)
#' relativeLogExpression(bcb_small)
NULL



# matrix =======================================================================
.relativeLogExpression.matrix <-  # nolint
    function(object) {
        t(t(object) / colMedians(object))
    }



#' @rdname relativeLogExpression
#' @export
setMethod(
    f = "relativeLogExpression",
    signature = signature("matrix"),
    definition = .relativeLogExpression.matrix
)



# SummarizedExperiment =========================================================
.relativeLogExpression.SE <-  # nolint
    function(object) {
        relativeLogExpression(counts(object))
    }



#' @rdname relativeLogExpression
#' @export
setMethod(
    f = "relativeLogExpression",
    signature = signature("SummarizedExperiment"),
    definition = .relativeLogExpression.SE
)
