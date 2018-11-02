# FIXME Move generics to basejump.generics



#' @rdname alphaSummary
#' @export
setGeneric(
    name = "alphaSummary",
    def = function(object, ...) {
        standardGeneric("alphaSummary")
    }
)



#' @rdname contrastName
#' @export
setGeneric(
    name = "contrastName",
    def = function(object, ...) {
        standardGeneric("contrastName")
    }
)



#' @rdname DESeqResultsTables
#' @export
setGeneric(
    name = "DESeqResultsTables",
    def = function(object, ...) {
        standardGeneric("DESeqResultsTables")
    }
)



#' @rdname plot5Prime3PrimeBias
#' @export
setGeneric(
    name = "plot5Prime3PrimeBias",
    def = function(object, ...) {
        standardGeneric("plot5Prime3PrimeBias")
    }
)



#' @rdname plotDEGHeatmap
#' @export
setGeneric(
    name = "plotDEGHeatmap",
    def = function(object, ...) {
        standardGeneric("plotDEGHeatmap")
    }
)



#' @rdname plotDEGPCA
#' @export
setGeneric(
    name = "plotDEGPCA",
    def = function(object, ...) {
        standardGeneric("plotDEGPCA")
    }
)



#' @rdname plotExonicMappingRate
#' @export
setGeneric(
    name = "plotExonicMappingRate",
    def = function(object, ...) {
        standardGeneric("plotExonicMappingRate")
    }
)



#' @rdname plotGeneSaturation
#' @export
setGeneric(
    name = "plotGeneSaturation",
    def = function(object, ...) {
        standardGeneric("plotGeneSaturation")
    }
)



#' @rdname plotIntronicMappingRate
#' @export
setGeneric(
    name = "plotIntronicMappingRate",
    def = function(object, ...) {
        standardGeneric("plotIntronicMappingRate")
    }
)



#' @rdname plotMappedReads
#' @export
setGeneric(
    name = "plotMappedReads",
    def = function(object, ...) {
        standardGeneric("plotMappedReads")
    }
)



#' @rdname plotMappingRate
#' @export
setGeneric(
    name = "plotMappingRate",
    def = function(object, ...) {
        standardGeneric("plotMappingRate")
    }
)



#' @rdname plotMeanSD
#' @export
setGeneric(
    name = "plotMeanSD",
    def = function(object, ...) {
        standardGeneric("plotMeanSD")
    }
)



#' @rdname plotRRNAMappingRate
#' @export
setGeneric(
    name = "plotRRNAMappingRate",
    def = function(object, ...) {
        standardGeneric("plotRRNAMappingRate")
    }
)



#' @rdname plotTotalReads
#' @export
setGeneric(
    name = "plotTotalReads",
    def = function(object, ...) {
        standardGeneric("plotTotalReads")
    }
)



#' @rdname plotVolcano
#' @export
setGeneric(
    name = "plotVolcano",
    def = function(object, ...) {
        standardGeneric("plotVolcano")
    }
)



# Don't use "rle" because that's for run length encoding.
#' @rdname relativeLogExpression
#' @export
setGeneric(
    name = "relativeLogExpression",
    def = function(object, ...) {
        standardGeneric("relativeLogExpression")
    }
)



#' @rdname tmm
#' @export
setGeneric(
    name = "tmm",
    def = function(object, ...) {
        standardGeneric("tmm")
    }
)



#' @rdname topTables
#' @export
setGeneric(
    name = "topTables",
    def = function(object, ...) {
        standardGeneric("topTables")
    }
)
