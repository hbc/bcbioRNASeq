## FIXME Need to define in AcidGenerics, so we can also use in DESeqAnalysis.
#' @rdname coerce
#' @export
setGeneric(
    name = "as.DESeqDataSet",
    def = function(x, ...) {
        standardGeneric("as.DESeqDataSet")
    }
)



## FIXME Need to define in AcidGenerics, so we can also use in DESeqAnalysis.
#' @rdname coerce
#' @export
setGeneric(
    name = "as.DESeqTransform",
    def = function(x, ...) {
        standardGeneric("as.DESeqTransform")
    }
)



## FIXME Need to define in AcidGenerics, so we can also use in DESeqAnalysis.
#' @rdname coerce
#' @export
setGeneric(
    name = "as.DGEList",
    def = function(x, ...) {
        standardGeneric("as.DGEList")
    }
)



#' @rdname plot5Prime3PrimeBias
#' @name plot5Prime3PrimeBias
#' @importFrom AcidGenerics plot5Prime3PrimeBias
#' @usage plot5Prime3PrimeBias(object, ...)
#' @export
NULL



#' @rdname plotPseudoVsAlignedCounts
#' @name plotPseudoVsAlignedCounts
#' @importFrom AcidGenerics plotPseudoVsAlignedCounts
#' @usage plotPseudoVsAlignedCounts(object, ...)
#' @export
NULL



#' @rdname slotAlignedCounts
#' @name slotAlignedCounts
#' @importFrom AcidGenerics slotAlignedCounts
#' @usage slotAlignedCounts(object, ...)
#' @export
NULL
