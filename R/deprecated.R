# nocov start
# nolint start



#' @name defunct
#' @inherit basejump::defunct
#' @keywords internal
NULL

#' @name deprecated
#' @inherit basejump::deprecated
#' @keywords internal
NULL



# v0.0.25 ======================================================================
#' @rdname deprecated
#' @export
plotGeneDetectionSaturation <- function(...) {
    .Deprecated("plotGeneSaturation")
    plotGeneSaturation(...)
}

#' @rdname deprecated
#' @export
plotDispersion <- function(...) {
    .Deprecated("plotDispEsts")
    plotDispEsts(...)
}



# v0.0.27 ======================================================================
#' @rdname deprecated
#' @export
loadRNASeqRun <- function(...) {
    .Deprecated("loadRNASeq")
    loadRNASeq(...)
}



# v0.1.2 =======================================================================
#' @rdname defunct
#' @export
plotGeneHeatmap <- function(...) {
    .Defunct("plotHeatmap")
}



# v0.1.3 =======================================================================
#' @rdname defunct
#' @export
txi <- function(...) {
    .Defunct("assays")
}



# v0.2.0 =======================================================================
# annotable deprecation for SummarizedExperiment added to bcbioBase v0.2.0

#' @rdname deprecated
#' @export
meltLog10 <- function(object, ...) {
    .Defunct("meltCounts")
}

#' @rdname deprecated
#' @export
setMethod(
    "design",
    signature("bcbioRNASeq"),
    function(object, ...) {
        .Defunct(msg = "Object no longer supports design formula")
    }
)

#' @rdname deprecated
#' @importFrom BiocGenerics design<-
#' @export
setMethod(
    "design<-",
    signature(
        object = "bcbioRNASeq",
        value = "formula"
    ),
    function(object, ..., value) {
        .Defunct(msg = "Object no longer supports design formula")
    }
)



# v0.2.2 =======================================================================
#' @rdname deprecated
#' @export
loadRNASeq <- function(...) {
    .Deprecated("bcbioRNASeq")
    bcbioRNASeq(...)
}



# v0.3.0 =======================================================================
#' @rdname defunct
#' @export
plotCountDensity <- function(...) {
    .Defunct("plotCountsPerGene(object, geom = \"density\")")
}

#' @rdname defunct
#' @export
resultsTables <- function(...) {
    .Defunct("DESeqResultsTables")
}



# nolint end
# nocov end
