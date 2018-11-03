# TODO Consider moving this to basejump.experiment.



#' @name relativeLogExpression
#' @inherit basejump.generics::relativeLogExpression
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams params
#'
#' @examples
#' data(bcb)
#' relativeLogExpression(bcb)
NULL



#' @importFrom basejump.generics relativeLogExpression
#' @aliases NULL
#' @export
basejump.generics::relativeLogExpression



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
