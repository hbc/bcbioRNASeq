#' @rdname coerce
#' @export
setGeneric(
    name = "as.DESeqDataSet",
    def = function(x, ...) {
        standardGeneric("as.DESeqDataSet")
    }
)



#' @rdname coerce
#' @export
setGeneric(
    name = "as.DESeqTransform",
    def = function(x, ...) {
        standardGeneric("as.DESeqTransform")
    }
)



#' @rdname coerce
#' @export
setGeneric(
    name = "as.DGEList",
    def = function(x, ...) {
        standardGeneric("as.DGEList")
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
