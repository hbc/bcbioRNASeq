#' Sample Metrics
#'
#' @rdname metrics
#' @name metrics
#'
#' @importFrom basejump metrics
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
#' @importFrom dplyr left_join
#' @importFrom magrittr set_rownames
#' @export
setMethod(
    "metrics",
    signature("bcbioRNASeq"),
    function(object) {
        suppressMessages(left_join(
            as.data.frame(colData(object)),
            as.data.frame(metadata(object)[["metrics"]])
        )) %>%
            set_rownames(.[["sampleID"]])
    })
