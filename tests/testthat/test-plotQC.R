context("plotQC")

test_that("plotQC", {
    x <- plotQC(object)
    expect_s3_class(x, "ggplot")
})

with_parameters_test_that(
    "Interesting groups support", {
        expect_is(fun, "function")
        x <- fun(object)
        expect_s3_class(x, "ggplot")
        x <- fun(object, interestingGroups = "sampleName")
        expect_s3_class(x, "ggplot")
    },
    fun = list(
        ## plotCountsPerGene,
        ## plotGenesDetected,
        plot5Prime3PrimeBias,
        plotCountsPerFeature,
        plotExonicMappingRate,
        plotGeneSaturation,
        plotFeaturesDetected,
        plotIntronicMappingRate,
        plotMappedReads,
        plotMappingRate,
        plotPCA,
        plotRRNAMappingRate,
        plotTotalReads
    )
)
