#' @rdname plotCountsCorrelation
#' @export
setGeneric(
    name = "plotCountsCorrelation",
    def = function(x, y, ...) {
        standardGeneric("plotCountsCorrelation")
    }
)



#' @rdname plotPseudoVsAlignedCounts
#' @export
setGeneric(
    name = "plotPseudoVsAlignedCounts",
    def = function(object, ...) {
        standardGeneric("plotPseudoVsAlignedCounts")
    }
)



#' @rdname slotAlignedCounts
#' @export
setGeneric(
    name = "slotAlignedCounts",
    def = function(object, ...) {
        standardGeneric("slotAlignedCounts")
    }
)
