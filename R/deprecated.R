## nocov start



#' @name deprecated
#' @inherit AcidRoxygen::deprecated description examples return seealso title
#' @inheritParams AcidRoxygen::params
#' @keywords internal
NULL



## v0.3.17 =====================================================================

#' @export
#' @rdname plotCountsPerFeature
plotCountsPerGene <- function(object, ...) {
    assert(.isGeneLevel(object))
    plotCountsPerFeature(object, ...)
}

#' @export
#' @rdname deprecated
plotGenesDetected <- function(object, ...) {
    assert(.isGeneLevel(object))
    plotFeaturesDetected(object, ...)
}



## F1000v2 manuscript ==========================================================

`plotDEGHeatmap,deprecated` <- # nolint
    function(object, results, counts, ...) {
        assert(
            is(results, "DESeqResults"),
            is(counts, "DESeqTransform")
        )
        plotDEGHeatmap(
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
    f = "plotDEGHeatmap",
    signature = signature(object = "missing"),
    definition = `plotDEGHeatmap,deprecated`
)

#' @export
#' @rdname deprecated
plotMeanAverage <- function(...) {
    ## > .Deprecated("plotMA")
    assert(requireNamespace("BiocGenerics", quietly = TRUE))
    BiocGenerics::plotMA(...)
}

#' @export
#' @rdname deprecated
prepareRNASeqTemplate <- function(...) {
    .Defunct("prepareTemplate")
}

`topTables,DataFrameList` <- # nolint
    function(object, ...) {
        markdownTables(
            object = object,
            ...
        )
    }

#' @export
#' @rdname deprecated
setMethod(
    f = "topTables",
    signature = signature(object = "DataFrameList"),
    definition = `topTables,DataFrameList`
)

#' @export
#' @rdname deprecated
writeCounts <-
    function(
        ...,
        dir = getOption(x = "acid.export.dir", default = getwd())
    ) {
    ## > .Deprecated("export")
    objects <- list(...)
    names(objects) <- dots(..., character = TRUE)
    Map(
        object = objects,
        con = file.path(dir, paste0(names(objects), ".csv")),
        f = export
    )
}



## nocov end
