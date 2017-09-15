#' Sample Metrics
#'
#' @rdname metrics
#' @name metrics
#'
#' @return [data.frame].
#'
#' @examples
#' data(bcb)
#' metrics(bcb) %>% glimpse
NULL



# Methods ====
#' @rdname metrics
#' @export
setMethod("metrics", "bcbioRNADataSet", function(object) {
    left_join(
        as.data.frame(colData(object)),
        as.data.frame(metadata(object)[["metrics"]]),
        by = c("sampleID", "sampleName")
    )
})
