## nocov start
## nolint start



#' @name deprecated
#' @inherit AcidRoxygen::deprecated description examples return seealso title
#' @inheritParams AcidRoxygen::params
#' @keywords internal
NULL



## v0.3.17 =====================================================================
#' @rdname plotCountsPerFeature
#' @export
plotCountsPerGene <- function(object, ...) {
    assert(.isGeneLevel(object))
    plotCountsPerFeature(object, ...)
}

#' @rdname deprecated
#' @export
plotGenesDetected <- function(object, ...) {
    assert(.isGeneLevel(object))
    plotFeaturesDetected(object, ...)
}



## nolint end
## nocov end
