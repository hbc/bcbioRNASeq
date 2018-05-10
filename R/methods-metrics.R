#' Sample Metrics
#'
#' @name metrics
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase metrics
#'
#' @inheritParams general
#'
#' @return `data.frame`.
#'
#' @examples
#' metrics(bcb_small) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("bcbioRNASeq"),
    function(object, interestingGroups) {
        validObject(object)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        } else {
            interestingGroups(object) <- interestingGroups
        }
        data <- as.data.frame(colData(object))
        data <- uniteInterestingGroups(data, interestingGroups)
        data
    }
)
