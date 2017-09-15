#' S4 Generics
#'
#' @rdname AllGenerics
#' @name AllGenerics
#' @keywords internal
#'
#' @param object Object.
#' @param x Object.
#' @param i An integer or numeric scalar.
#' @param withDimnames A `logical`, indicating whether dimnames should be
#'   applied to extracted assay elements.
#' @param ... *Additional arguments (for the S4 generic definition).*
#'
#' @return No value.
NULL



#' Quality Control Plots
#'
#' @rdname qcPlots
#' @name qcPlots
#' @keywords internal
#'
#' @param counts Object containing a count matrix.
#' @param flip Flip X and Y axes.
#' @param interestingGroup Category to use to group samples (color and shape).
#'   If unset, this is automatically determined by the metadata set inside the
#'   [bcbioRNADataSet].
#' @param minCounts Numeric value for filtering the counts matrix before
#'   plotting.
#' @param normalized Count normalization method. See [counts()] documentation
#'   for more information.
#' @param passLimit Threshold to plot pass color marker.
#' @param warnLimit Threshold to plot warning color marker.
#'
#' @return [ggplot].
NULL



#' @rdname alphaSummary
#' @family Differential Expression Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("alphaSummary", function(object, ...) {
    standardGeneric("alphaSummary")
})



#' @rdname meltLog10
#' @inheritParams AllGenerics
#' @export
setGeneric("meltLog10", function(object, ...) {
    standardGeneric("meltLog10")
})



#' @rdname plot53Bias
#' @inheritParams AllGenerics
#' @export
setGeneric("plot53Bias", function(object, ...) {
    standardGeneric("plot53Bias")
})



#' @rdname plotCorrelationHeatmap
#' @inheritParams AllGenerics
#' @export
setGeneric("plotCorrelationHeatmap", function(object, ...) {
    standardGeneric("plotCorrelationHeatmap")
})



#' @rdname plotCountDensity
#' @inheritParams AllGenerics
#' @export
setGeneric("plotCountDensity", function(object, ...) {
    standardGeneric("plotCountDensity")
})



#' @rdname plotCountsPerGene
#' @inheritParams AllGenerics
#' @export
setGeneric("plotCountsPerGene", function(object, ...) {
    standardGeneric("plotCountsPerGene")
})



#' @rdname plotDEGHeatmap
#' @inheritParams AllGenerics
#' @export
setGeneric("plotDEGHeatmap", function(object, counts, ...) {
    standardGeneric("plotDEGHeatmap")
})



#' @rdname plotExonicMappingRate
#' @inheritParams AllGenerics
#' @export
setGeneric("plotExonicMappingRate", function(object, ...) {
    standardGeneric("plotExonicMappingRate")
})



#' @rdname plotGenderMarkers
#' @inheritParams AllGenerics
#' @export
setGeneric("plotGenderMarkers", function(object, ...) {
    standardGeneric("plotGenderMarkers")
})



#' @rdname plotGene
#' @inheritParams AllGenerics
#' @export
setGeneric("plotGene", function(object, ...) {
    standardGeneric("plotGene")
})



#' @rdname plotGeneHeatmap
#' @inheritParams AllGenerics
#' @export
setGeneric("plotGeneHeatmap", function(object, ...) {
    standardGeneric("plotGeneHeatmap")
})



#' @rdname plotGeneSaturation
#' @inheritParams AllGenerics
#' @export
setGeneric(
    "plotGeneSaturation",
    function(object, counts, ...) {
        standardGeneric("plotGeneSaturation")
    })



#' @rdname plotGenesDetected
#' @inheritParams AllGenerics
#' @export
setGeneric("plotGenesDetected", function(object, counts, ...) {
    standardGeneric("plotGenesDetected")
})



#' @rdname plotIntronicMappingRate
#' @inheritParams AllGenerics
#' @export
setGeneric("plotIntronicMappingRate", function(object, ...) {
    standardGeneric("plotIntronicMappingRate")
})



#' @rdname plotMappedReads
#' @inheritParams AllGenerics
#' @export
setGeneric("plotMappedReads", function(object, ...) {
    standardGeneric("plotMappedReads")
})



#' @rdname plotMappingRate
#' @inheritParams AllGenerics
#' @export
setGeneric("plotMappingRate", function(object, ...) {
    standardGeneric("plotMappingRate")
})



#' @rdname plotMeanSD
#' @inheritParams AllGenerics
#' @export
setGeneric("plotMeanSD", function(object, ...) {
    standardGeneric("plotMeanSD")
})



#' @rdname plotPCACovariates
#' @inheritParams AllGenerics
#' @export
setGeneric("plotPCACovariates", function(object, ...) {
    standardGeneric("plotPCACovariates")
})



#' @rdname plotRRNAMappingRate
#' @inheritParams AllGenerics
#' @export
setGeneric("plotRRNAMappingRate", function(object, ...) {
    standardGeneric("plotRRNAMappingRate")
})



#' @rdname plotTotalReads
#' @inheritParams AllGenerics
#' @export
setGeneric("plotTotalReads", function(object, ...) {
    standardGeneric("plotTotalReads")
})



#' @rdname plotVolcano
#' @inheritParams AllGenerics
#' @export
setGeneric("plotVolcano", function(object, ...) {
    standardGeneric("plotVolcano")
})



#' @rdname prepareRNASeqTemplate
#' @inheritParams AllGenerics
#' @export
setGeneric("prepareRNASeqTemplate", function(object, ...) {
    standardGeneric("prepareRNASeqTemplate")
})



#' @rdname resultsTables
#' @inheritParams AllGenerics
#' @export
setGeneric("resultsTables", function(object, ...) {
    standardGeneric("resultsTables")
})



#' @rdname tmm
#' @inheritParams AllGenerics
#' @export
setGeneric("tmm", function(object) {
    standardGeneric("tmm")
})



#' @rdname topTables
#' @inheritParams AllGenerics
#' @export
setGeneric("topTables", function(object, ...) {
    standardGeneric("topTables")
})



#' @rdname tpm
#' @inheritParams AllGenerics
#' @export
setGeneric("tpm", function(object) {
    standardGeneric("tpm")
})



#' @rdname txi
#' @inheritParams AllGenerics
#' @export
setGeneric("txi", function(object) {
    standardGeneric("txi")
})
