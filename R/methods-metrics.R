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
        sampleData <- as.data.frame(sampleData(object))
        metrics <- metadata(object)[["metrics"]]
        assert_are_identical(rownames(sampleData), rownames(metrics))
        assert_are_disjoint_sets(colnames(sampleData), colnames(metrics))
        cbind(sampleData, metrics) %>%
            as.data.frame() %>%
            rownames_to_column() %>%
            mutate_if(is.character, as.factor) %>%
            column_to_rownames()
    }
)
