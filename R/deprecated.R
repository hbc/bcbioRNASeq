#' Deprecated Functions
#'
#' @rdname deprecated
#' @name deprecated
#' @keywords internal
#' @inheritParams AllGenerics
#' @return No value.
NULL



# v0.0.25 ======================================================================
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



# v0.0.27 ======================================================================
#' @rdname deprecated
#' @export
loadRNASeqRun <- function(...) {
    .Deprecated("loadRNASeq")
    loadRNASeq(...)
}



# v0.1.2 =======================================================================
#' @rdname deprecated
#' @export
plotGeneHeatmap <- function(...) {
    .Deprecated("plotHeatmap")
    plotHeatmap(...)
}



# v0.1.3 =======================================================================
#' @rdname deprecated
#' @export
txi <- function(object) {
    .Deprecated("bcbio(object, \"tximport\")")
    bcbio(object, "tximport")
}
