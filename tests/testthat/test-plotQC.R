test_that("plotQC", {
    x <- plotQC(object)
    expect_s3_class(x, "ggplot")
})

test_that("Fast mode", {
    x <- plotQC(bcb_fast)
    expect_s3_class(x, "ggplot")
})

test_that("Interesting groups support", {
    for (fun in list(
        ## plotCountsPerGene
        ## plotGenesDetected
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
    )) {
        expect_is(fun, "function")
        x <- fun(object)
        expect_s3_class(x, "ggplot")
        x <- fun(object, interestingGroups = "sampleName")
        expect_s3_class(x, "ggplot")
    }
})
