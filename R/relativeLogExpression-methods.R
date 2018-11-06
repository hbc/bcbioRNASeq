# TODO Consider moving this to basejump.



#' @name relativeLogExpression
#' @inherit basejump::relativeLogExpression
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams params
#'
#' @examples
#' data(bcb)
#' relativeLogExpression(bcb)
NULL



#' @importFrom basejump relativeLogExpression
#' @aliases NULL
#' @export
basejump::relativeLogExpression



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
