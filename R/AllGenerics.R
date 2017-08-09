#' S4 Generics
#'
#' @rdname AllGenerics
#' @name AllGenerics
#' @keywords internal
#'
#' @param object Object.
#' @param i An integer or numeric scalar.
#' @param withDimnames A `logical`, indicating whether dimnames should be
#'   applied to extracted assay elements.
#' @param x Object.
#' @param ... Additional arguments.
#'
#' @return No value.
NULL



#' Quality Control Plots
#'
#' @rdname qcPlots
#' @name qcPlots
#' @keywords internal
#'
#' @param counts Counts matrix.
#' @param flip Flip X and Y axes.
#' @param interestingGroup *Optional*. Category to use to group samples (color
#'   and shape). If unset, this is automatically determined by the metadata set
#'   inside the [bcbioRNADataSet].
#' @param minCounts Numeric value for filtering the counts matrix before
#'   plotting.
#' @param normalized Count normalization method. See [counts()] documentation
#'   for more information.
#' @param passLimit Threshold to plot pass color marker.
#' @param warnLimit Threshold to plot warning color marker.
#'
#' @return [ggplot].
NULL



#' @rdname aggregateReplicates
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("aggregateReplicates", function(object, ...) {
    standardGeneric("aggregateReplicates")
})



#' @rdname alphaSummary
#' @author Michael Steinbaugh, Lorena Patano
#' @inheritParams AllGenerics
#' @export
setGeneric("alphaSummary", function(object, ...) {
    standardGeneric("alphaSummary")
})



#' @rdname bcbio
#' @author Lorena Pantano, Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("bcbio", function(object, ...) {
    standardGeneric("bcbio")
})



#' @rdname bcbio
#' @export
setGeneric("bcbio<-", function(object, ..., value) {
    standardGeneric("bcbio<-")
})



#' @rdname download
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("download", function(object) {
    standardGeneric("download")
})



#' @rdname loadRNASeqRun
#' @author Michael Steinbaugh, Lorena Pantano
#' @inheritParams AllGenerics
#' @export
setGeneric("loadRNASeqRun", function(object, ...) {
    standardGeneric("loadRNASeqRun")
})



#' @rdname meltLog10
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("meltLog10", function(object, ...) {
    standardGeneric("meltLog10")
})



#' @rdname metadataTable
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("metadataTable", function(object, ...) {
    standardGeneric("metadataTable")
})



#' @rdname metrics
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("metrics", function(object) {
    standardGeneric("metrics")
})



#' @rdname plotCorrelationHeatmap
#' @family Heatmaps
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("plotCorrelationHeatmap", function(object, ...) {
    standardGeneric("plotCorrelationHeatmap")
})



#' @rdname plotCountDensity
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inheritParams AllGenerics
#' @inherit qcPlots
#' @export
setGeneric("plotCountDensity", function(object, ...) {
    standardGeneric("plotCountDensity")
})



#' @rdname plotCountsPerGene
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inheritParams AllGenerics
#' @inherit qcPlots
#' @export
setGeneric("plotCountsPerGene", function(object, ...) {
    standardGeneric("plotCountsPerGene")
})



#' @rdname plotDEGHeatmap
#' @family Heatmaps
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("plotDEGHeatmap", function(object, counts, ...) {
    standardGeneric("plotDEGHeatmap")
})



#' @rdname plotDispersion
#' @family DESeq2 Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("plotDispersion", function(object) {
    standardGeneric("plotDispersion")
})



#' @rdname plotExonicMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inheritParams AllGenerics
#' @inherit qcPlots
#' @export
setGeneric("plotExonicMappingRate", function(object, ...) {
    standardGeneric("plotExonicMappingRate")
})



#' @rdname plotGenderMarkers
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("plotGenderMarkers", function(object, ...) {
    standardGeneric("plotGenderMarkers")
})



#' @rdname plotGene
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("plotGene", function(object, ...) {
    standardGeneric("plotGene")
})



#' @rdname plotGeneHeatmap
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("plotGeneHeatmap", function(object, ...) {
    standardGeneric("plotGeneHeatmap")
})



#' @rdname plotGeneSaturation
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inheritParams AllGenerics
#' @inherit qcPlots
#' @export
setGeneric(
    "plotGeneSaturation",
    function(object, counts, ...) {
        standardGeneric("plotGeneSaturation")
    })



#' @rdname plotGenesDetected
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inheritParams AllGenerics
#' @inherit qcPlots
#' @export
setGeneric("plotGenesDetected", function(object, counts, ...) {
    standardGeneric("plotGenesDetected")
})



#' @rdname plotIntronicMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inheritParams AllGenerics
#' @inherit qcPlots
#' @export
setGeneric("plotIntronicMappingRate", function(object, ...) {
    standardGeneric("plotIntronicMappingRate")
})



#' @rdname plotMappedReads
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inheritParams AllGenerics
#' @inherit qcPlots
#' @export
setGeneric("plotMappedReads", function(object, ...) {
    standardGeneric("plotMappedReads")
})



#' @rdname plotMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inheritParams AllGenerics
#' @inherit qcPlots
#' @export
setGeneric("plotMappingRate", function(object, ...) {
    standardGeneric("plotMappingRate")
})



#' @rdname plotMeanSD
#' @family DESeq2 Utilities
#' @author Michael Steinbaugh, Lorena Patano
#' @inheritParams AllGenerics
#' @export
setGeneric("plotMeanSD", function(object, ...) {
    standardGeneric("plotMeanSD")
})



#' @rdname plotPCACovariates
#' @author Lorena Pantano, Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("plotPCACovariates", function(object, ...) {
    standardGeneric("plotPCACovariates")
})



#' @rdname plotRRNAMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inheritParams AllGenerics
#' @inherit qcPlots
#' @export
setGeneric("plotRRNAMappingRate", function(object, ...) {
    standardGeneric("plotRRNAMappingRate")
})



#' @rdname plotTotalReads
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inheritParams AllGenerics
#' @inherit qcPlots
#' @export
setGeneric("plotTotalReads", function(object, ...) {
    standardGeneric("plotTotalReads")
})



#' @rdname plotVolcano
#' @family Differential Expression Plots
#' @author John Hutchinson, Michael Steinbaugh, Lorena Pantano
#' @inheritParams AllGenerics
#' @export
setGeneric("plotVolcano", function(object, ...) {
    standardGeneric("plotVolcano")
})



#' @rdname plot53Bias
#' @family Quality Control Plots
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @inherit qcPlots
#' @export
setGeneric("plot53Bias", function(object, ...) {
    standardGeneric("plot53Bias")
})



#' @rdname resultsTables
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("resultsTables", function(object, ...) {
    standardGeneric("resultsTables")
})



#' @rdname sampleDirs
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("sampleDirs", function(object) {
    standardGeneric("sampleDirs")
})



#' @rdname tmm
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("tmm", function(object) {
    standardGeneric("tmm")
})



#' @rdname topTables
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("topTables", function(object, ...) {
    standardGeneric("topTables")
})



#' @rdname tpm
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("tpm", function(object) {
    standardGeneric("tpm")
})



#' @rdname txi
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("txi", function(object) {
    standardGeneric("txi")
})
