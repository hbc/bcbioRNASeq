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
        # Stop on fast-rnaseq pipline detection
        if (!"totalReads" %in% colnames(colData(object))) {
            # Parse the YAML metadata or log file here instead?
            stop(paste(
                "Fast mode detected.",
                "Metrics were not calculated."
            ), call. = FALSE)
        }
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
