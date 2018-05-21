context("Quality Control Functions")

skip <- bcb_small
assays(skip)[["rlog"]] <- NULL
assays(skip)[["vst"]] <- NULL

skipWarning <- paste(
    "rlog not present in assays.",
    "Calculating log2 TMM counts instead."
)



# Metrics ======================================================================
test_that("Quality Control Metrics Plots", {
    fxns <- c(
        "plot5Prime3PrimeBias",
        "plotCountDensity",
        "plotCountsPerGene",
        "plotExonicMappingRate",
        "plotGeneSaturation",
        "plotGenesDetected",
        "plotIntronicMappingRate",
        "plotMappedReads",
        "plotMappingRate",
        "plotMeanSD",
        "plotPCA",
        "plotRRNAMappingRate",
        "plotTotalReads"
    )
    invisible(lapply(fxns, function(fxn) {
        fxn <- get(fxn, inherits = TRUE)
        expect_is(fxn, "nonstandardGenericFunction")
        p <- fxn(bcb_small)
        expect_is(p, "ggplot")
    }))
})



# plotCorrelationHeatmap =======================================================
test_that("plotCorrelationHeatmap", {
    # Pearson (default)
    p <- plotCorrelationHeatmap(bcb_small)
    expect_identical(names(p), pheatmapNames)
    # Spearman
    p <- plotCorrelationHeatmap(bcb_small, method = "spearman")
    expect_identical(names(p), pheatmapNames)
    # Bad method
    expect_error(
        plotCorrelationHeatmap(bcb_small, method = "XXX"),
        "'arg' should be one of"
    )
})

test_that("plotCorrelationHeatmap : transformationLimit", {
    expect_warning(
        plotCorrelationHeatmap(skip, normalized = "rlog"),
        skipWarning
    )
    p <- suppressWarnings(plotCorrelationHeatmap(skip, normalized = "rlog"))
    expect_identical(names(p), pheatmapNames)
})



# plotPCA ======================================================================
test_that("plotPCA : transformationLimit", {
    expect_warning(
        plotPCA(skip, normalized = "rlog"),
        skipWarning
    )
    p <- suppressWarnings(plotPCA(skip, normalized = "rlog"))
    expect_is(p, "ggplot")
})



# plotPCACovariates ============================================================
test_that("plotPCACovariates", {
    # BioC 3.6 version of DEGreport returns warnings.
    # Test against GitHub version >= 1.15
    p <- plotPCACovariates(bcb_small)
    expect_is(p, "list")
    expect_identical(
        names(p),
        c("significantCovars",
          "plot",
          "corMatrix",
          "pcsMatrix",
          "scatterPlot",
          "effectsSignificantCovars")
    )
})

test_that("plotPCACovariates : Significant covars", {
    p <- plotPCACovariates(
        bcb_small,
        metrics = c("exonicRate", "intronicRate")
    )
    # Don't expect these to be significant with the example dataset
    expect_identical(
        as.character(p[["significantCovars"]]),
        character()
    )
})

test_that("plotPCACovariates : Skip DESeq2 transforms", {
    expect_warning(
        plotPCACovariates(skip, normalized = "rlog"),
        skipWarning
    )
    p <- suppressWarnings(plotPCA(skip, normalized = "rlog"))
    expect_is(p, "ggplot")
})



# plotDispEsts =================================================================
test_that("plotDispEsts", {
    p <- plotDispEsts(bcb_small)
    expect_is(p, "list")
    expect_identical(
        names(p),
        c("rect", "text")
    )
})
