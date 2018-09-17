context("Quality Control Functions")

skip <- bcb_small
assays(skip)[["rlog"]] <- NULL
assays(skip)[["vst"]] <- NULL

skipWarning <- paste(
    "rlog not present in assays.",
    "Calculating log2 TMM counts instead."
)



# QC plots with interesting groups =============================================
with_parameters_test_that(
    "Plots supporting interesting groups", {
        expect_is(fun, "function")
        object <- fun(object = bcb_small)
        expect_is(object, "ggplot")
        object <- fun(object = bcb_small, interestingGroups = "sampleName")
        expect_is(object, "ggplot")
    },
    fun = list(
        plot5Prime3PrimeBias,
        plotCountDensity,
        plotCountsPerGene,
        plotExonicMappingRate,
        plotGeneSaturation,
        plotGenesDetected,
        plotIntronicMappingRate,
        plotMappedReads,
        plotMappingRate,
        plotPCA,
        plotRRNAMappingRate,
        plotTotalCounts,
        plotTotalReads
    )
)



# plotCorrelationHeatmap =======================================================
with_parameters_test_that(
    "plotCorrelationHeatmap", {
        # Pearson
        expect_is(
            object = plotCorrelationHeatmap(object, method = "pearson"),
            class = "pheatmap"
        )
        # Spearman
        expect_is(
            object = plotCorrelationHeatmap(object, method = "spearman"),
            class = "pheatmap"
        )
        # Bad method
        expect_error(
            object = plotCorrelationHeatmap(object, method = "XXX"),
            regexp = "'arg' should be one of"
        )
    },
    object = list(
        bcbioRNASeq = bcb_small,
        DESeqDataSet = dds_small,
        DESeqResults = vst_small
    )
)

test_that("plotCorrelationHeatmap : Skipped DESeq transform", {
    expect_warning(
        object = plotCorrelationHeatmap(skip, normalized = "rlog"),
        regexp = skipWarning
    )
    expect_is(
        object = suppressWarnings(
            plotCorrelationHeatmap(skip, normalized = "rlog")
        ),
        class = "pheatmap"
    )
})



# plotCountsPerGene ============================================================
test_that("plotCountsPerGene", {
    expect_is(
        object = plotCountsPerGene(
            object = bcb_small,
            normalized = "vst",
            title = NULL,
            interestingGroups = "sampleName"
        ),
        class = "ggplot"
    )
})



# plotCountDensity =============================================================
test_that("plotCountDensity", {
    # solid style
    object <- plotCountDensity(bcb_small, normalized = "tmm", style = "solid")
    expect_is(object, "ggplot")

    # vst
    # Set title = NULL to check disabling of subtitle
    object <- plotCountDensity(
        object = bcb_small,
        normalized = "vst",
        interestingGroups = "sampleName",
        title = NULL
    )
    expect_is(object, "ggplot")
})



# plotGeneSaturation ===========================================================
test_that("plotGeneSaturation", {
    object <- plotGeneSaturation(
        object = bcb_small,
        trendline = TRUE,
        label = TRUE
    )
    expect_is(object, "ggplot")
})



# plotMeanSD ===================================================================
test_that("plotMeanSD : DESeqDataSet", {
    object <- plotMeanSD(dds_small)
    expect_is(object, "ggplot")
})

test_that("plotMeanSD : bcbioRNASeq : No stashed DESeq transforms", {
    x <- bcb_small
    assays(x)[["rlog"]] <- NULL
    assays(x)[["vst"]] <- NULL
    object <- plotMeanSD(x)
    expect_is(object, "ggplot")
})



# plotPCA ======================================================================
test_that("plotPCA : Label", {
    object <- plotPCA(bcb_small, label = FALSE)
    expect_is(object, "ggplot")
    object <- plotPCA(bcb_small, label = TRUE)
    expect_is(object, "ggplot")
})

test_that("plotPCA : data.frame", {
    object <- plotPCA(bcb_small, return = "data.frame")
    expect_is(object, "data.frame")
})

test_that("plotPCA : transformationLimit", {
    expect_warning(
        plotPCA(skip, normalized = "rlog"),
        skipWarning
    )
    object <- suppressWarnings(plotPCA(skip, normalized = "rlog"))
    expect_is(object, "ggplot")
})



# plotPCACovariates ============================================================
test_that("plotPCACovariates", {
    # BioC 3.6 version of DEGreport returns warnings.
    # Test against GitHub version >= 1.15
    object <- plotPCACovariates(bcb_small)
    expect_is(object, "list")
    expect_identical(
        names(object),
        c("significantCovars",
          "plot",
          "corMatrix",
          "pcsMatrix",
          "scatterPlot",
          "effectsSignificantCovars")
    )
})

test_that("plotPCACovariates : Metrics", {
    object <- plotPCACovariates(
        object = bcb_small,
        metrics = c("exonicRate", "intronicRate")
    )
    # Don't expect these to be significant with the example dataset
    expect_identical(
        as.character(object[["significantCovars"]]),
        character()
    )

    # If metrics = FALSE, we require at least 2 interesting groups
    expect_error(
        object <- plotPCACovariates(bcb_small, metrics = FALSE),
        "`plotPCACovariates\\(\\)` requires >= 2 metrics"
    )
})

test_that("plotPCACovariates : Skip DESeq2 transforms", {
    expect_warning(
        plotPCACovariates(skip, normalized = "rlog"),
        skipWarning
    )
    object <- suppressWarnings(plotPCA(skip, normalized = "rlog"))
    expect_is(object, "ggplot")
})



# plotDispEsts =================================================================
test_that("plotDispEsts", {
    object <- plotDispEsts(bcb_small)
    expect_is(object, "list")
    expect_identical(
        names(object),
        c("rect", "text")
    )
})
