#' @rdname alphaSummary
#' @export
setGeneric(
    "alphaSummary",
    function(object, ...) {
        standardGeneric("alphaSummary")
    }
)



#' @rdname contrastName
#' @export
setGeneric(
    "contrastName",
    function(object, ...) {
        standardGeneric("contrastName")
    }
)



#' @rdname plot5Prime3PrimeBias
#' @export
setGeneric(
    "plot5Prime3PrimeBias",
    function(object, ...) {
        standardGeneric("plot5Prime3PrimeBias")
    }
)



#' @rdname plotCountDensity
#' @export
setGeneric(
    "plotCountDensity",
    function(object, ...) {
        standardGeneric("plotCountDensity")
    }
)



#' @rdname plotCountsPerGene
#' @export
setGeneric(
    "plotCountsPerGene",
    function(object, ...) {
        standardGeneric("plotCountsPerGene")
    }
)



#' @rdname plotDEGHeatmap
#' @export
setGeneric(
    "plotDEGHeatmap",
    function(results, counts, ...) {
        standardGeneric("plotDEGHeatmap")
    }
)



#' @rdname plotDEGPCA
#' @export
setGeneric(
    "plotDEGPCA",
    function(results, counts, ...) {
        standardGeneric("plotDEGPCA")
    }
)



#' @rdname plotExonicMappingRate
#' @export
setGeneric(
    "plotExonicMappingRate",
    function(object, ...) {
        standardGeneric("plotExonicMappingRate")
    }
)



#' @rdname plotGenderMarkers
#' @export
setGeneric(
    "plotGenderMarkers",
    function(object, ...) {
        standardGeneric("plotGenderMarkers")
    }
)



#' @rdname plotGeneSaturation
#' @export
setGeneric(
    "plotGeneSaturation",
    function(object, ...) {
        standardGeneric("plotGeneSaturation")
    }
)



#' @rdname plotGenesDetected
#' @export
setGeneric(
    "plotGenesDetected",
    function(object, ...) {
        standardGeneric("plotGenesDetected")
    }
)



#' @rdname plotIntronicMappingRate
#' @export
setGeneric(
    "plotIntronicMappingRate",
    function(object, ...) {
        standardGeneric("plotIntronicMappingRate")
    }
)



# Alternate version, so we don't conflict with `DESeq2::plotMA` methods
#' @rdname plotMeanAverage
#' @export
setGeneric(
    "plotMeanAverage",
    function(object, ...) {
        standardGeneric("plotMeanAverage")
    }
)



#' @rdname plotMappedReads
#' @export
setGeneric(
    "plotMappedReads",
    function(object, ...) {
        standardGeneric("plotMappedReads")
    }
)



#' @rdname plotMappingRate
#' @export
setGeneric(
    "plotMappingRate",
    function(object, ...) {
        standardGeneric("plotMappingRate")
    }
)



#' @rdname plotMeanSD
#' @export
setGeneric(
    "plotMeanSD",
    function(object, ...) {
        standardGeneric("plotMeanSD")
    }
)



#' @rdname plotPCACovariates
#' @export
setGeneric(
    "plotPCACovariates",
    function(object, ...) {
        standardGeneric("plotPCACovariates")
    }
)



#' @rdname plotRRNAMappingRate
#' @export
setGeneric(
    "plotRRNAMappingRate",
    function(object, ...) {
        standardGeneric("plotRRNAMappingRate")
    }
)



#' @rdname plotTotalReads
#' @export
setGeneric(
    "plotTotalReads",
    function(object, ...) {
        standardGeneric("plotTotalReads")
    }
)



#' @rdname plotVolcano
#' @export
setGeneric(
    "plotVolcano",
    function(object, ...) {
        standardGeneric("plotVolcano")
    }
)



#' @rdname resultsTables
#' @export
setGeneric(
    "resultsTables",
    function(results, counts, ...) {
        standardGeneric("resultsTables")
    }
)



#' @rdname tmm
#' @export
setGeneric(
    "tmm",
    function(object) {
        standardGeneric("tmm")
    }
)



#' @rdname topTables
#' @export
setGeneric(
    "topTables",
    function(object, ...) {
        standardGeneric("topTables")
    }
)



#' @rdname tpm
#' @export
setGeneric(
    "tpm",
    function(object) {
        standardGeneric("tpm")
    }
)
