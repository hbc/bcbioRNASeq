#' Sample Metrics
#'
#' @name metrics
#' @family Data Functions
#' @author Michael Steinbaugh
#' @importFrom basejump metrics
#' @export
#'
#' @inheritParams general
#'
#' @return `data.frame`.
#'
#' @examples
#' x <- metrics(bcb_small)
#' glimpse(x)
NULL



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("bcbioRNASeq"),
    function(object, interestingGroups = NULL) {
        validObject(object)
        if (!is.null(interestingGroups)) {
            interestingGroups(object) <- interestingGroups
        }
        data <- do.call(
            what = metrics,
            args = list(
                object = as(object, "SummarizedExperiment")
            )
        )
        # Returning as `data.frame` instead of `DataFrame`.
        # `metrics()` return is being used in tidyverse/ggplot scripts.
        # Consider returning `DataFrame` instead in a future update.
        as.data.frame(data)
    }
)
