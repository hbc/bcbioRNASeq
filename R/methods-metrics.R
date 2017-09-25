#' Sample Metrics
#'
#' @rdname metrics
#' @name metrics
#'
#' @inheritParams AllGenerics
#'
#' @return [data.frame].
#'
#' @examples
#' data(bcb)
#' metrics(bcb) %>% str
NULL



# Methods ====
#' @rdname metrics
#' @export
setMethod("metrics", "bcbioRNADataSet", function(object) {
    left_join(
        as.data.frame(colData(object)),
        as.data.frame(metadata(object)[["metrics"]]),
        by = c("sampleID", "sampleName")
    ) %>%
        set_rownames(.[["sampleID"]])
})
