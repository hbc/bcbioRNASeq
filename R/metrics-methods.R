#' @name metrics
#' @inherit bioverbs::metrics
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param ... Additional arguments.
#'
#' @return `data.frame`.
#'
#' @examples
#' metrics(bcb_small) %>% glimpse()
NULL



#' @rdname metrics
#' @name metrics
#' @importFrom bioverbs metrics
#' @usage metrics(object, ...)
#' @export
NULL



metrics.bcbioRNASeq <-  # nolint
    function(object, interestingGroups) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )

        # Stop on fast-rnaseq pipline detection.
        # Consider parsing the YAML metadata as an alternate approach.
        if (!"totalReads" %in% colnames(colData(object))) {
            # nocov start
            stop(paste(
                "Fast mode detected.",
                "Metrics were not calculated."
            ), call. = FALSE)
            # nocov end
        }

        data <- as.data.frame(colData(object))
        data <- uniteInterestingGroups(data, interestingGroups)

        assert_is_subset(
            x = c("sampleName", "interestingGroups"),
            y = colnames(data)
        )
        data
    }



#' @rdname metrics
#' @export
setMethod(
    f = "metrics",
    signature = signature("bcbioRNASeq"),
    definition = metrics.bcbioRNASeq
)
