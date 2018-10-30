context("Plots : Quality Control")

data(bcb, deseq, envir = environment())
dds_small <- as(deseq, "DESeqDataSet")
vst_small <- as(deseq, "DESeqTransform")



# QC plots with interesting groups =============================================
with_parameters_test_that(
    "Plots supporting interesting groups", {
        expect_is(fun, "function")
        object <- fun(object = bcb)
        expect_is(object, "ggplot")
        object <- fun(object = bcb, interestingGroups = "sampleName")
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
        bcbioRNASeq = bcb,
        DESeqDataSet = dds_small,
        DESeqResults = vst_small
    )
)



# plotCountsPerGene ============================================================
# FIXME Can we parameterize the functions that support normalized?
test_that("plotCountsPerGene", {
    expect_is(
        object = plotCountsPerGene(
            object = bcb,
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
        x <- plotCountsPerGene(bcb, geom = geom)
        expect_is(x, "ggplot")

        x <- plotCountsPerGene(
            object = bcb,
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
        object = bcb,
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
        bcbioRNASeq = bcb,
        DESeqDataSet = dds_small
    )
)



# plotPCA ======================================================================
test_that("plotPCA : Label", {
    object <- plotPCA(bcb, label = FALSE)
    expect_is(object, "ggplot")
    object <- plotPCA(bcb, label = TRUE)
    expect_is(object, "ggplot")
})

test_that("plotPCA : DataFrame", {
    object <- plotPCA(bcb, return = "DataFrame")
    expect_is(object, "DataFrame")
})



# plotDispEsts =================================================================
test_that("plotDispEsts", {
    object <- plotDispEsts(bcb)
    expect_is(object, "list")
    expect_identical(
        names(object),
        c("rect", "text")
    )
})



context("Plots : Gene Expression")

g2s <- Gene2Symbol(bcb)
geneIDs <- head(g2s[["geneID"]])
geneNames <- head(g2s[["geneName"]])



# plotGenderMarkers ============================================================
with_parameters_test_that(
    "plotGenderMarkers", {
        expect_is(
            object = plotGenderMarkers(object),
            class = "ggplot"
        )
    },
    object = list(
        bcbioRNASeq = bcb,
        DESeqDataSet = dds_small,
        DESeqTransform = vst_small
    )
)



# plotGene =====================================================================
with_parameters_test_that(
    "plotGene", {
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
},
    object = list(
        bcbioRNASeq = bcb,
        DESeqDataSet = dds_small,
        DESeqTransform = vst_small,
        DESeqAnalysis = deseq
    )
)



# plotHeatmap ==================================================================
test_that("plotHeatmap : bcbioRNASeq", {
    genes <- head(rownames(bcb), n = 100L)
    p <- plotHeatmap(bcb[genes, ])
    expect_is(p, "pheatmap")
})

test_that("plotHeatmap : DESeqDataSet", {
    genes <- head(rownames(dds_small), n = 20L)
    p <- plotHeatmap(dds_small[genes, ])
    expect_is(p, "pheatmap")
})
