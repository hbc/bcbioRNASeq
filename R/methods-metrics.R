#' Sample Metrics
#'
#' @rdname metrics
#' @author Michael Steinbaugh
#'
#' @return [tibble].
#' @export
#'
#' @examples
#' data(bcb)
#' metrics(bcb)
setMethod("metrics", "bcbioRNADataSet", function(object) {
    metrics <- metadata(object)[["metrics"]]
    if (is.null(metrics)) return(NULL)
    col_data <- colData(object) %>% as("tibble")
    left_join(col_data, metrics, by = meta_priority_cols)
})
