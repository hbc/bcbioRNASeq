#' Deprecated Functions
#'
#' @rdname deprecated
#' @name deprecated
#'
#' @return No value.
NULL



#' @rdname deprecated
#' @export
download <- function(...) {
    .Deprecated("prepareRNASeqTemplate")
    prepareRNASeqTemplate(...)
}
