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



#' @rdname plotPCA
#' @name plotPCA
#' @importFrom AcidGenerics plotPCA
#' @usage plotPCA(object, ...)
#' @export
NULL



#' @rdname plotPCACovariates
#' @name plotPCACovariates
#' @importFrom AcidGenerics plotPCACovariates
#' @usage plotPCACovariates(object, ...)
#' @export
NULL



#' @rdname plotPseudoVsAlignedCounts
#' @name plotPseudoVsAlignedCounts
#' @importFrom AcidGenerics plotPseudoVsAlignedCounts
#' @usage plotPseudoVsAlignedCounts(object, ...)
#' @export
NULL



#' @rdname plotRRNAMappingRate
#' @name plotRRNAMappingRate
#' @importFrom AcidGenerics plotRRNAMappingRate
#' @usage plotRRNAMappingRate(object, ...)
#' @export
NULL



#' @rdname plotTotalReads
#' @name plotTotalReads
#' @importFrom AcidGenerics plotTotalReads
#' @usage plotTotalReads(object, ...)
#' @export
NULL



#' @rdname relativeLogExpression
#' @name relativeLogExpression
#' @importFrom AcidGenerics relativeLogExpression
#' @usage relativeLogExpression(object, ...)
#' @export
NULL



## FIXME show



#' @rdname slotAlignedCounts
#' @name slotAlignedCounts
#' @importFrom AcidGenerics slotAlignedCounts
#' @usage slotAlignedCounts(object, ...)
#' @export
NULL



#' @rdname tmm
#' @name tmm
#' @importFrom AcidGenerics tmm
#' @usage tmm(object, ...)
#' @export
NULL



## FIXME updateObject
