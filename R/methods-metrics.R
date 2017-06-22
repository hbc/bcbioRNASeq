#' Sample metrics
#'
#' @rdname metrics
#' @docType methods
#'
#' @param object [bcbioRnaDataSet] object.
#'
#' @export
setMethod("metrics", "bcbioRnaDataSet", function(object) {
    metrics <- metadata(object)[["metrics"]]
    if (is.null(metrics)) return(NULL)
    cbind(colData(object), metrics) %>% as.data.frame
})
