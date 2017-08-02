#' Sample Metrics
#'
#' @rdname metrics
#' @author Michael Steinbaugh
#'
#' @return [data.frame].
#' @export
#'
#' @examples
#' data(bcb)
#' metrics(bcb)
setMethod("metrics", "bcbioRNADataSet", function(object) {
    metrics <- .uniqueMetrics(object)
    if (is.null(metrics)) return(NULL)
    meta <- .interestingColData(object)
    cbind(meta, metrics)
})
