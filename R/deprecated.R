## nocov start
## nolint start



#' @name defunct
#' @inherit basejump::defunct
#' @keywords internal
NULL



#' @name deprecated
#' @inherit basejump::deprecated
#' @keywords internal
NULL



## v0.2.2 =======================================================================
#' @rdname defunct
#' @export
loadRNASeq <- function(...) {
    .Defunct("bcbioRNASeq")
}



## v0.3.0 =======================================================================
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



## v0.3.16 ======================================================================
#' @rdname deprecated
#' @export
prepareRNASeqTemplate <- function(...) {
    .Deprecated("prepareTemplate(package = \"bcbioRNASeq\")")
    prepareTemplate(package = "bcbioRNASeq", ...)
}



## v0.3.17 ======================================================================
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



## v0.3.22 =====================================================================
#' @importFrom basejump aggregateReplicates
#' @export
basejump::aggregateReplicates

#' @rdname deprecated
#' @export
alphaSummary <- function(object, ...) {
    ## > .Deprecated("DESeqAnalysis::alphaSummary")
    ## This doesn't work unless we attach the package.
    require("DESeqAnalysis", quietly = TRUE)
    DESeqAnalysis::alphaSummary(object, ...)
}

#' @importFrom DESeqAnalysis contrastName
#' @export
DESeqAnalysis::contrastName

#' @rdname deprecated
#' @export
plotDEGHeatmap <- function(results, counts, ...) {
    ## > .Deprecated("DESeqAnalysis::plotDEGHeatmap")
    requireNamespace("DESeqAnalysis", quietly = TRUE)
    DESeqAnalysis::plotDEGHeatmap(
        object = results,
        counts = counts,
        ...
    )
}

#' @rdname deprecated
#' @export
plotDEGPCA <- function(results, counts, ...) {
    ## > .Deprecated("DESeqAnalysis::plotDEGPCA")
    requireNamespace("DESeqAnalysis", quietly = TRUE)
    DESeqAnalysis::plotDEGPCA(
        object = results,
        counts = counts,
        ...
    )
}

#' @importFrom DESeqAnalysis plotMA
#' @export
DESeqAnalysis::plotMA

#' @importFrom DESeqAnalysis plotMeanAverage
#' @export
DESeqAnalysis::plotMeanAverage

#' @importFrom DESeqAnalysis plotVolcano
#' @export
DESeqAnalysis::plotVolcano

## FIXME Add support for counts in DESeqAnalysis...
## FIXME This doesn't support Dropbox any more. Consider updating the method
## to work with previous bcbioRNASeq support.

#' @rdname deprecated
#' @export
resultsTables <- function(results, counts, ...) {
    ## > .Deprecated("DESeqAnalysis::plotDEGPCA")
    resultsTables(
        object = results,
        counts = counts,
        ...
    )
}

#' @importFrom DESeqAnalysis topTables
#' @export
DESeqAnalysis::topTables



## nolint end
## nocov end
