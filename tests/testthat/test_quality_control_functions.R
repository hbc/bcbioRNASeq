context("Quality Control Functions")

bcb_skip <- bcb_small
assays(bcb_skip)[["rlog"]] <- NULL
assays(bcb_skip)[["vst"]] <- NULL

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
        plotCountsPerGene,
        plotExonicMappingRate,
        plotGeneSaturation,
        plotGenesDetected,
        plotIntronicMappingRate,
        plotMappedReads,
        plotMappingRate,
        plotPCA,
        plotRRNAMappingRate,
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
        object = plotCorrelationHeatmap(bcb_skip, normalized = "rlog"),
        regexp = skipWarning
    )
    expect_is(
        object = suppressWarnings(
            plotCorrelationHeatmap(bcb_skip, normalized = "rlog")
        ),
        class = "pheatmap"
    )
})



# plotCountsPerGene ============================================================
# FIXME Can we parameterize the functions that support normalized?
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



# plotCountsPerGene ============================================================
with_parameters_test_that(
    "plotCountsPerGene", {
        x <- plotCountsPerGene(bcb_small, geom = geom)
        expect_is(x, "ggplot")

        x <- plotCountsPerGene(
            object = bcb_small,
            normalized = "vst",
            interestingGroups = "sampleName",
            title = NULL
        )
        expect_is(x, "ggplot")
    },
    geom = methodFormals(
        f = "plotCountsPerGene",
        signature = "bcbioRNASeq"
    ) %>%
        .[["geom"]] %>%
        eval()
)



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
with_parameters_test_that(
    "plotMeanSD", {
        x <- plotMeanSD(object)
        expect_is(x, "ggplot")
    },
    object = list(
        bcbioRNASeq1 = bcb_small,
        bcbioRNASeq2 = bcb_skip,
        DESeqDataSet = dds_small
    )
)



# plotPCA ======================================================================
test_that("plotPCA : Label", {
    object <- plotPCA(bcb_small, label = FALSE)
    expect_is(object, "ggplot")
    object <- plotPCA(bcb_small, label = TRUE)
    expect_is(object, "ggplot")
})

test_that("plotPCA : DataFrame", {
    object <- plotPCA(bcb_small, return = "DataFrame")
    expect_is(object, "DataFrame")
})

test_that("plotPCA : Skipped DESeq2 transforms", {
    expect_warning(
        plotPCA(bcb_skip, normalized = "rlog"),
        skipWarning
    )
    object <- suppressWarnings(plotPCA(bcb_skip, normalized = "rlog"))
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
