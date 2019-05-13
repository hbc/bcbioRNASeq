context("Quality Control Functions")

skip <- bcb_small
assays(skip)[["rlog"]] <- NULL
assays(skip)[["vst"]] <- NULL

skipWarning <- paste(
    "rlog not present in assays.",
    "Calculating log2 TMM counts instead."
)



# Loop across QC plots =========================================================
test_that("Plots supporting interesting groups", {
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
        "plotPCA",
        "plotRRNAMappingRate",
        "plotTotalReads"
    )
    invisible(lapply(fxns, function(fxn) {
        fxn <- get(fxn, inherits = TRUE)
        expect_is(fxn, "nonstandardGenericFunction")
        p <- fxn(object = bcb_small)
        expect_is(p, "ggplot")
        p <- fxn(object = bcb_small, interestingGroups = "sampleName")
        expect_is(p, "ggplot")
    }))
})



# plotCorrelationHeatmap =======================================================
test_that("plotCorrelationHeatmap : bcbioRNASeq", {
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

test_that("plotCorrelationHeatmap : DESeqTransform", {
    p <- plotCorrelationHeatmap(vst_small)
    expect_identical(names(p), pheatmapNames)
})

test_that("plotCorrelationHeatmap : transformationLimit", {
    expect_warning(
        plotCorrelationHeatmap(skip, normalized = "rlog"),
        skipWarning
    )
    p <- suppressWarnings(plotCorrelationHeatmap(skip, normalized = "rlog"))
    expect_identical(names(p), pheatmapNames)
})



# plotCountsPerGene ============================================================
test_that("plotCountsPerGene", {
    p <- plotCountsPerGene(
        object = bcb_small,
        normalized = "vst",
        title = NULL,
        interestingGroups = "sampleName"
    )
    expect_is(p, "ggplot")
})



# plotCountDensity =============================================================
test_that("plotCountDensity", {
    # solid style
    p <- plotCountDensity(bcb_small, normalized = "tmm", style = "solid")
    expect_is(p, "ggplot")

    # vst
    # Set title = NULL to check disabling of subtitle
    p <- plotCountDensity(
        object = bcb_small,
        normalized = "vst",
        interestingGroups = "sampleName",
        title = NULL
    )
    expect_is(p, "ggplot")
})



# plotGeneSaturation ===========================================================
test_that("plotGeneSaturation", {
    p <- plotGeneSaturation(
        object = bcb_small,
        trendline = TRUE,
        label = TRUE
    )
    expect_is(p, "ggplot")
})



# plotMeanSD ===================================================================
test_that("plotMeanSD : DESeqDataSet", {
    p <- plotMeanSD(dds_small)
    expect_is(p, "ggplot")
})

test_that("plotMeanSD : bcbioRNASeq : No stashed DESeq transforms", {
    x <- bcb_small
    assays(x)[["rlog"]] <- NULL
    assays(x)[["vst"]] <- NULL
    p <- plotMeanSD(x)
    expect_is(p, "ggplot")
})



# plotPCA ======================================================================
test_that("plotPCA : Label", {
    p <- plotPCA(bcb_small, label = FALSE)
    expect_is(p, "ggplot")
    p <- plotPCA(bcb_small, label = TRUE)
    expect_is(p, "ggplot")
})

test_that("plotPCA : data.frame", {
    p <- plotPCA(bcb_small, return = "data.frame")
    expect_is(p, "data.frame")
})

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
    p <- plotPCACovariates(bcb_small)
    expect_is(p, "list")
    expect_identical(
        names(p),
        c(
            "plot",
            "corMatrix",
            "pcsMatrix",
            "scatterPlot",
            "significants"
        )
    )
})

test_that("plotPCACovariates : Metrics", {
    p <- plotPCACovariates(
        object = bcb_small,
        metrics = c("exonicRate", "intronicRate")
    )
    # Don't expect these to be significant with the example dataset.
    expect_identical(
        as.character(p[["significantCovars"]]),
        character()
    )

    # If metrics = FALSE, we require at least 2 interesting groups.
    expect_error(
        object = plotPCACovariates(bcb_small, metrics = FALSE),
        regexp = "2 metrics"
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
