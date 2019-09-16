## nocov start
## nolint start



#' @name defunct
#' @inherit acidroxygen::defunct description examples return seealso title
#' @inheritParams acidroxygen::params
#' @keywords internal
NULL



#' @name deprecated
#' @inherit acidroxygen::deprecated description examples return seealso title
#' @inheritParams acidroxygen::params
#' @keywords internal
NULL



## v0.2.2 =======================================================================
#' @rdname defunct
#' @export
loadRNASeq <- function(...) {
    .Defunct("bcbioRNASeq")
}



## v0.3.16 ======================================================================
#' @rdname defunct
#' @export
prepareRNASeqTemplate <- function(...) {
    .Defunct("prepareTemplate(package = \"bcbioRNASeq\")")
}



## v0.3.17 ======================================================================
#' @rdname deprecated
#' @export
plotCountsPerGene <- function(object, ...) {
    .Deprecated("plotCountsPerFeature")
    assert(.isGeneLevel(object))
    plotCountsPerFeature(object, ...)
}

#' @rdname deprecated
#' @export
plotGenesDetected <- function(object, ...) {
    .Deprecated("plotFeaturesDetected")
    assert(.isGeneLevel(object))
    plotFeaturesDetected(object, ...)
}



## nolint end
## nocov end
