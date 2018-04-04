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
    function(object) {
        validObject(object)
        metrics <- metadata(object)[["metrics"]]
        metadata <- sampleData(object)
        assert_are_identical(rownames(metrics), rownames(metadata))

        if (length(intersect(colnames(metrics), legacyMetricsCols))) {
            warn(paste(
                paste(
                    "Metrics slot contains legacy columns:",
                    toString(intersect(colnames(metrics), legacyMetricsCols))
                ),
                bcbioBase::updateMessage,
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
