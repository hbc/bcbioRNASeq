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
    .Defunct("plotCountsPerFeature(object, geom = \"density\")")
}

#' @rdname defunct
#' @export
resultsTables <- function(...) {
    .Defunct(msg = paste(
        "resultsTables approach has been reworked in DESeqAnalysis package.",
        "https://deseqanalysis.acidgenomics.com/",
        sep = "\n"
    ))
}



# v0.3.12 ======================================================================
# SummarizedExperiment method now works seamlessly with bcbioRNASeq object, so
# no need to redefine a custom method here.
#' @importFrom basejump sampleData
#' @export
basejump::sampleData



# v0.3.16 ======================================================================
#' @rdname deprecated
#' @export
prepareRNASeqTemplate <- function(...) {
    .Deprecated("prepareTemplate(package = \"bcbioRNASeq\")")
    prepareTemplate(package = "bcbioRNASeq", ...)
}



# v0.3.17 ======================================================================
#' @rdname deprecated
#' @export
plotCountsPerGene <- function(object, title = "Counts per gene", ...) {
    assert(.isGeneLevel(object))
    do.call(
        what = plotCountsPerFeature,
        args = matchArgsToDoCall()
    )
}

#' @rdname deprecated
#' @export
plotGenesDetected <- function(object, ...) {
    assert(.isGeneLevel(object))
    plotFeaturesDetected(
        object = object,
        countsAxisLabel = "genes",
        title = "Genes detected",
        ...
    )
}



# nolint end
# nocov end
