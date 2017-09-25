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
#' @export
setGeneric("alphaSummary", function(object, ...) {
    standardGeneric("alphaSummary")
})



#' @rdname meltLog10
#' @export
setGeneric("meltLog10", function(object, ...) {
    standardGeneric("meltLog10")
})



#' @rdname plot53Bias
#' @export
setGeneric("plot53Bias", function(object, ...) {
    standardGeneric("plot53Bias")
})



#' @rdname plotCorrelationHeatmap
#' @export
setGeneric("plotCorrelationHeatmap", function(object, ...) {
    standardGeneric("plotCorrelationHeatmap")
})



#' @rdname plotCountDensity
#' @export
setGeneric("plotCountDensity", function(object, ...) {
    standardGeneric("plotCountDensity")
})



#' @rdname plotCountsPerGene
#' @export
setGeneric("plotCountsPerGene", function(object, ...) {
    standardGeneric("plotCountsPerGene")
})



#' @rdname plotDEGHeatmap
#' @export
setGeneric("plotDEGHeatmap", function(object, counts, ...) {
    standardGeneric("plotDEGHeatmap")
})



#' @rdname plotExonicMappingRate
#' @export
setGeneric("plotExonicMappingRate", function(object, ...) {
    standardGeneric("plotExonicMappingRate")
})



#' @rdname plotGenderMarkers
#' @export
setGeneric("plotGenderMarkers", function(object, ...) {
    standardGeneric("plotGenderMarkers")
})



#' @rdname plotGeneHeatmap
#' @export
setGeneric("plotGeneHeatmap", function(object, ...) {
    standardGeneric("plotGeneHeatmap")
})



#' @rdname plotGeneSaturation
#' @export
setGeneric(
    "plotGeneSaturation",
    function(object, counts, ...) {
        standardGeneric("plotGeneSaturation")
    })



#' @rdname plotGenesDetected
#' @export
setGeneric("plotGenesDetected", function(object, counts, ...) {
    standardGeneric("plotGenesDetected")
})



#' @rdname plotIntronicMappingRate
#' @export
setGeneric("plotIntronicMappingRate", function(object, ...) {
    standardGeneric("plotIntronicMappingRate")
})



#' @rdname plotMappedReads
#' @export
setGeneric("plotMappedReads", function(object, ...) {
    standardGeneric("plotMappedReads")
})



#' @rdname plotMappingRate
#' @export
setGeneric("plotMappingRate", function(object, ...) {
    standardGeneric("plotMappingRate")
})



#' @rdname plotMeanSD
#' @export
setGeneric("plotMeanSD", function(object, ...) {
    standardGeneric("plotMeanSD")
})



#' @rdname plotPCACovariates
#' @export
setGeneric("plotPCACovariates", function(object, ...) {
    standardGeneric("plotPCACovariates")
})



#' @rdname plotRRNAMappingRate
#' @export
setGeneric("plotRRNAMappingRate", function(object, ...) {
    standardGeneric("plotRRNAMappingRate")
})



#' @rdname plotTotalReads
#' @export
setGeneric("plotTotalReads", function(object, ...) {
    standardGeneric("plotTotalReads")
})



#' @rdname plotVolcano
#' @export
setGeneric("plotVolcano", function(object, ...) {
    standardGeneric("plotVolcano")
})



#' @rdname prepareRNASeqTemplate
#' @export
setGeneric("prepareRNASeqTemplate", function(object, ...) {
    standardGeneric("prepareRNASeqTemplate")
})



#' @rdname resultsTables
#' @export
setGeneric("resultsTables", function(object, ...) {
    standardGeneric("resultsTables")
})



#' @rdname tmm
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
#' @export
setGeneric("tpm", function(object) {
    standardGeneric("tpm")
})



#' @rdname txi
#' @export
setGeneric("txi", function(object) {
    standardGeneric("txi")
})
