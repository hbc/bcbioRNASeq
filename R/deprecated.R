## nocov start
## nolint start



#' @name defunct
#' @inherit AcidRoxygen::defunct description examples return seealso title
#' @inheritParams AcidRoxygen::params
#' @keywords internal
NULL



#' @name deprecated
#' @inherit AcidRoxygen::deprecated description examples return seealso title
#' @inheritParams AcidRoxygen::params
#' @keywords internal
NULL



## v0.3.16 =====================================================================
#' @rdname defunct
#' @export
prepareRNASeqTemplate <- function(...) {
    .Defunct("prepareTemplate(package = \"bcbioRNASeq\")")
}



## v0.3.17 =====================================================================
#' @rdname plotCountsPerFeature
#' @export
plotCountsPerGene <- function(object, ...) {
    assert(.isGeneLevel(object))
    plotCountsPerFeature(object, ...)
}

#' @rdname deprecated
#' @export
plotGenesDetected <- function(object, ...) {
    assert(.isGeneLevel(object))
    plotFeaturesDetected(object, ...)
}



## v0.3.29 =====================================================================
## Check if DESeqAnalysis package is installed, otherwise inform user about
## deprecation of DESeq2 methods inside bcbioRNASeq package.
.requireDESeqAnalysis <- function(name) {
    if (!requireNamespace("DESeqAnalysis")) {
        stop(
            "'", name, "()' has migrated to DESeqAnalysis package.\n",
            "Run 'BiocManager::install(\"acidgenomics/DESeqAnalysis\")' ",
            "to install."
        )
    }
}

#' @rdname deprecated
#' @export
alphaSummary <- function(...) {
    .requireDESeqAnalysis(name = "alphaSummary")
    DESeqAnalysis::alphaSummary(...)
}

#' @rdname deprecated
#' @export
contrastName <- function(...) {
    .requireDESeqAnalysis(name = "contrastName")
    DESeqAnalysis::contrastName(...)
}

#' @rdname deprecated
#' @export
plotDEGHeatmap <- function(...) {
    .requireDESeqAnalysis(name = "plotDEGHeatmap")
    DESeqAnalysis::plotDEGHeatmap(...)
}

#' @rdname deprecated
#' @export
plotDEGPCA <- function(...) {
    .requireDESeqAnalysis(name = "plotDEGPCA")
    DESeqAnalysis::plotDEGPCA(...)
}

#' @rdname deprecated
#' @export
plotMA <- function(...) {
    .requireDESeqAnalysis(name = "plotMA")
    DESeqAnalysis::plotMA(...)
}

#' @rdname deprecated
#' @export
plotMeanAverage <- function(...) {
    .requireDESeqAnalysis(name = "plotMeanAverage")
    DESeqAnalysis::plotMeanAverage(...)
}

#' @rdname deprecated
#' @export
plotVolcano <- function(...) {
    .requireDESeqAnalysis(name = "plotVolcano")
    DESeqAnalysis::plotVolcano(...)
}

#' @rdname deprecated
#' @export
resultsTables <- function(...) {
    .requireDESeqAnalysis(name = "resultsTables")
    DESeqAnalysis::resultsTables(...)
}

#' @rdname deprecated
#' @export
topTables <- function(...) {
    .requireDESeqAnalysis(name = "topTables")
    DESeqAnalysis::topTables(...)
}



## nolint end
## nocov end
