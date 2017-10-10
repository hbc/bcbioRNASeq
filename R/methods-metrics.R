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
#' metrics(bcb) %>%
#'     str()
NULL



# Methods ====
#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("bcbioRNASeqANY"),
    function(object) {
        suppressMessages(left_join(
            as.data.frame(colData(object)),
            as.data.frame(metadata(object)[["metrics"]])
        )) %>%
            set_rownames(.[["sampleID"]])
    })
