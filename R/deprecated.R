## nocov start



#' @name deprecated
#' @inherit AcidRoxygen::deprecated description examples return seealso title
#' @inheritParams AcidRoxygen::params
#' @keywords internal
NULL



#' @export
#' @rdname deprecated
camel <- function(...) {
    .Deprecated("syntactic::camelCase")
    assert(requireNamespaces("syntactic"))
    syntactic::camelCase(...)
}



#' @export
#' @rdname deprecated
plotCountsPerGene <- function(object, ...) {
    .Deprecated("plotCountsPerFeature")
    assert(.isGeneLevel(object))
    plotCountsPerFeature(object, ...)
}



#' @export
#' @rdname deprecated
plotDEGHeatmap <- function(...) {
    .Deprecated("DESeqAnalysis::plotDegHeatmap")
    assert(requireNamespaces("DESeqAnalysis"))
    DESeqAnalysis::plotDegHeatmap(...)
}



`plotDegHeatmap,deprecated` <- # nolint
    function(object, results, counts, ...) {
        assert(
            is(results, "DESeqResults"),
            is(counts, "DESeqTransform")
        )
        plotDegHeatmap(
            object = results,
            DESeqTransform = counts,
            ...
        )
    }

#' @export
#' @rdname deprecated
#' @param `results` `DESeqResults.`
#' @param `counts` `DESeqTransform`.
setMethod(
    f = "plotDegHeatmap",
    signature = signature(object = "missing"),
    definition = `plotDegHeatmap,deprecated`
)



#' @export
#' @rdname deprecated
plotGeneSaturation <- function(object, ...) {
    .Deprecated("plotFeatureSaturation")
    assert(.isGeneLevel(object))
    plotFeatureSaturation(object, ...)
}



#' @export
#' @rdname deprecated
plotGenesDetected <- function(object, ...) {
    .Deprecated("plotFeaturesDetected")
    assert(.isGeneLevel(object))
    plotFeaturesDetected(object, ...)
}



#' @export
#' @rdname deprecated
plotMA <- function(...) {
    .Deprecated("DESeqAnalysis::plotMa")
    assert(requireNamespaces("DESeqAnalysis"))
    DESeqAnalysis::plotMa(...)
}



#' @export
#' @rdname deprecated
plotMeanAverage <- function(...) {
    .Deprecated("DESeqAnalysis::plotMa")
    assert(requireNamespaces("DESeqAnalysis"))
    DESeqAnalysis::plotMa(...)
}



#' @export
#' @rdname deprecated
plotMeanSD <- function(...) {
    .Deprecated("plotMeanSd")
    plotMeanSd(...)
}



#' @export
#' @rdname deprecated
plotPCA <- function(...) {
    .Deprecated("plotPca")
    plotPca(...)
}



#' @export
#' @rdname deprecated
plotPCACovariates <- function(...) {
    .Deprecated("plotPcaCovariates")
    plotPcaCovariates(...)
}



#' @export
#' @rdname deprecated
plotQC <- function(...) {
    .Deprecated("plotQc")
    plotQc(...)
}



#' @export
#' @rdname deprecated
plotRRNAMappingRate <- function(...) {
    .Deprecated("plotRrnaMappingRate")
    plotRrnaMappingRate(...)
}



#' @export
#' @rdname deprecated
prepareRNASeqTemplate <- function(...) {
    .Deprecated("AcidMarkdown::prepareTemplate")
    assert(requireNamespaces("AcidMarkdown"))
    AcidMarkdown::prepareTemplate(...)
}



#' @export
#' @rdname deprecated
topTables <- function(...) {
    .Deprecated("DESeqAnalysis::markdownTables")
    assert(requireNamespaces("DESeqAnalysis"))
    DESeqAnalysis::markdownTables(...)
}



#' @export
#' @rdname deprecated
writeCounts <-
    function(...,
             dir = getOption(x = "acid.export.dir", default = getwd())) {
        .Deprecated("export")
        objects <- list(...)
        names(objects) <- dots(..., character = TRUE)
        out <- Map(
            object = objects,
            con = file.path(dir, paste0(names(objects), ".csv")),
            f = export
        )
        invisible(out)
    }



## nocov end
