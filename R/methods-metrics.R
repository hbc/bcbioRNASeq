#' Sample Metrics
#'
#' @rdname metrics
#' @name metrics
#'
#' @return [data.frame].
#' @export
#'
#' @examples
#' data(bcb)
#' metrics(bcb)
NULL



# Constructors ====
.metrics <- function(object) {
    metrics <- .uniqueMetrics(object)
    if (is.null(metrics)) return(NULL)
    meta <- .interestingColData(object)
    cbind(meta, metrics)
}



# Methods ====
#' @rdname metrics
#' @export
setMethod("metrics", "bcbioRNADataSet", .metrics)
