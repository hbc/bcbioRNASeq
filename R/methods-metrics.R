#' Sample Metrics
#'
#' @name metrics
#'
#' @importFrom bcbioBase metrics
#'
#' @inheritParams general
#'
#' @return `data.frame`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' metrics(bcb) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname metrics
#' @importFrom dplyr left_join mutate_if
#' @export
setMethod(
    "metrics",
    signature("bcbioRNASeq"),
    function(object) {
        validObject(object)
        metrics <- metadata(object)[["metrics"]]
        metadata <- colData(object, return = "data.frame")
        assert_are_identical(rownames(metrics), rownames(metadata))

        if (length(intersect(colnames(metrics), legacyMetricsCols))) {
            warn(paste(
                paste(
                    "Metrics slot contains legacy columns:",
                    toString(intersect(colnames(metrics), legacyMetricsCols))
                ),
                updateMsg,
                sep = "\n"
            ))
            metrics <- metrics %>%
                .[, setdiff(colnames(.), legacyMetricsCols), drop = FALSE]
        }
        assert_are_disjoint_sets(colnames(metrics), colnames(metadata))

        cbind(metadata, metrics) %>%
            rownames_to_column() %>%
            mutate_if(is.character, as.factor) %>%
            column_to_rownames()
    }
)
