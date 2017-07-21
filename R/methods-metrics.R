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
    metrics <- .unique_metrics(object)
    if (is.null(metrics)) return(NULL)
    meta <- .interesting_col_data(object)
    cbind(meta, metrics)
})
