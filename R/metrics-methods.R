#' Sample Metrics
#'
#' @name metrics
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase metrics
#' @export
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
        interestingGroups <- .prepareInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )

        # Stop on fast-rnaseq pipline detection
        if (!"totalReads" %in% colnames(colData(object))) {
            # Parse the YAML metadata or log file here instead?
            stop(paste(
                "Fast mode detected.",
                "Metrics were not calculated."
            ), call. = FALSE)
        }

        data <- as.data.frame(colData(object))
        data <- uniteInterestingGroups(data, interestingGroups)

        assert_is_subset(
            x = c("sampleName", "interestingGroups"),
            y = colnames(data)
        )
        data
    }
)
