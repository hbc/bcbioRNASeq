#' Sample metrics
#'
#' @rdname metrics
#' @docType methods
#'
#' @param object [bcbioRNADataSet] object.
#'
#' @export
setMethod("metrics", "bcbioRNADataSet", function(object) {
    metrics <- metadata(object)[["metrics"]]
    if (is.null(metrics)) return(NULL)
    cbind(colData(object), metrics) %>% as.data.frame
})
