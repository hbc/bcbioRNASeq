# FIXME Parameterize the functions that support `normalized` argument?

context("Plots : Quality Control")

data(bcb, envir = environment())
object <- bcb



# QC plots with interesting groups =============================================
with_parameters_test_that(
    "Plots supporting interesting groups", {
        expect_is(fun, "function")
        x <- fun(object)
        expect_s3_class(x, "ggplot")
        x <- fun(object, interestingGroups = "sampleName")
        expect_s3_class(x, "ggplot")
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
test_that("plotCorrelationHeatmap", {
    # Pearson.
    expect_s3_class(
        object = plotCorrelationHeatmap(object, method = "pearson"),
        class = "pheatmap"
    )
    # Spearman.
    expect_s3_class(
        object = plotCorrelationHeatmap(object, method = "spearman"),
        class = "pheatmap"
    )
    # Bad method.
    expect_error(
        object = plotCorrelationHeatmap(object, method = "XXX"),
        regexp = "'arg' should be one of"
    )
})



# plotCountsPerGene ============================================================
test_that("plotCountsPerGene", {
    expect_s3_class(
        object = plotCountsPerGene(
            object = object,
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
        x <- plotCountsPerGene(object, geom = geom)
        expect_s3_class(x, "ggplot")
        x <- plotCountsPerGene(
            object = object,
            normalized = "vst",
            interestingGroups = "sampleName",
            title = NULL
        )
        expect_s3_class(x, "ggplot")
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
        object = object,
        trendline = TRUE,
        label = TRUE
    )
    expect_s3_class(object, "ggplot")
})



# plotMeanSD ===================================================================
test_that("plotMeanSD", {
    x <- plotMeanSD(object)
    expect_s3_class(x, "ggplot")
})



# plotPCA ======================================================================
test_that("plotPCA : Label", {
    object <- plotPCA(bcb, label = FALSE)
    expect_s3_class(object, "ggplot")
    object <- plotPCA(bcb, label = TRUE)
    expect_s3_class(object, "ggplot")
})

test_that("plotPCA : DataFrame", {
    object <- plotPCA(bcb, return = "DataFrame")
    expect_s4_class(object, "DataFrame")
})



# plotDispEsts =================================================================
test_that("plotDispEsts", {
    object <- plotDispEsts(bcb)
    expect_type(object, "list")
    expect_identical(
        object = names(object),
        expected = c("rect", "text")
    )
})



context("Plots : Gene Expression")

g2s <- Gene2Symbol(bcb)
geneIDs <- head(g2s[["geneID"]])
geneNames <- head(g2s[["geneName"]])



# plotGenderMarkers ============================================================
test_that("plotGenderMarkers", {
    expect_s3_class(
        object = plotGenderMarkers(object),
        class = "ggplot"
    )
})



# plotGene =====================================================================
test_that("plotGene", {
    # Facet wrapped.
    p <- plotGene(
        object = object,
        genes = geneNames,
        interestingGroups = "sampleName",
        style = "facet"
    )
    expect_s3_class(p, "ggplot")

    # Wide format.
    p <- plotGene(
        object = bcb,
        genes = geneNames,
        style = "wide"
    )
    expect_s3_class(p, "ggplot")
})



# plotHeatmap ==================================================================
test_that("plotHeatmap", {
    genes <- head(rownames(bcb), n = 100L)
    p <- plotHeatmap(bcb[genes, ])
    expect_s3_class(p, "pheatmap")
})
