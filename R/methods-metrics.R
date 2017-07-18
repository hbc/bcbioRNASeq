#' Sample metrics
#'
#' @rdname metrics
#' @docType methods
#'
#' @param object [bcbioRNADataSet] object.
#'
#' @return [tibble].
#'
#' @export
setMethod("metrics", "bcbioRNADataSet", function(object) {
    metrics <- metadata(object)[["metrics"]]
    if (is.null(metrics)) return(NULL)
    col_data <- colData(object)
    left_join(col_data, metrics, by = meta_priority_cols) %>% as("tibble")
})
