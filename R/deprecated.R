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



## v0.3.0 =======================================================================
#' @rdname defunct
#' @export
plotCountDensity <- function(...) {
    .Defunct("plotCountsPerFeature(object, geom = \"density\")")
}



## v0.3.16 ======================================================================
#' @rdname deprecated
#' @export
prepareRNASeqTemplate <- function(...) {
    .Deprecated("prepareTemplate(package = \"bcbioRNASeq\")")
    prepareTemplate(package = "bcbioRNASeq", ...)
}



## v0.3.17 ======================================================================
#' @rdname deprecated
#' @export
plotCountsPerGene <- function(object, title = "Counts per gene", ...) {
    assert(.isGeneLevel(object))
    do.call(
        what = plotCountsPerFeature,
        args = matchArgsToDoCall()
    )
}

#' @rdname deprecated
#' @export
plotGenesDetected <- function(object, ...) {
    assert(.isGeneLevel(object))
    plotFeaturesDetected(
        object = object,
        countsAxisLabel = "genes",
        title = "Genes detected",
        ...
    )
}



## nolint end
## nocov end
