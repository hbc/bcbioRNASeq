#' Deprecated Functions
#'
#' @rdname deprecated
#' @name deprecated
#' @keywords internal
#' @inheritParams AllGenerics
#' @return No value.
NULL



#' @rdname deprecated
#' @export
download <- function(...) {
    .Deprecated("prepareRNASeqTemplate")
    prepareRNASeqTemplate(...)
}
