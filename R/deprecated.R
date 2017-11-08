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
download <- function() {
    .Deprecated("prepareRNASeqTemplate")
}

#' @rdname deprecated
#' @export
plotGeneDetectionSaturation <- function(...) {
    .Deprecated("plotGeneSaturation")
    plotGeneSaturation(...)
}

#' @rdname deprecated
#' @export
plotDispersion <- function(...) {
    .Deprecated("plotDispEsts")
    plotDispEsts(...)
}



# 0.0.27 ====
#' @rdname deprecated
#' @export
loadRNASeqRun <- function(...) {
    .Deprecated("loadRNASeq")
    loadRNASeq(...)
}



# 0.1.2 ====
#' @rdname deprecated
#' @export
plotGeneHeatmap <- function(...) {
    .Deprecated("plotHeatmap")
    plotHeatmap(...)
}
