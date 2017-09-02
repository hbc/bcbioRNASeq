#' Deprecated Functions
#'
#' @rdname deprecated
#' @name deprecated
#' @keywords internal
#' @inheritParams AllGenerics
#' @return No value.
NULL



# 0.0.25 ====
#' @rdname deprecated
#' @export
download <- function(...) {
    .Deprecated("prepareRNASeqTemplate")
    prepareRNASeqTemplate(...)
}

#' @rdname deprecated
#' @export
plotGeneDetectionSaturation <- function(...) {
    .Deprecated("plotGeneSaturation")
    plotGeneSaturation(...)
}
