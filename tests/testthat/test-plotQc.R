test_that("plotQc", {
    x <- plotQc(object)
    expect_s3_class(x, "ggplot")
})

test_that("Fast mode", {
    x <- plotQc(bcb_fast)
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
        plotPca,
        plotRrnaMappingRate,
        plotTotalReads
    )) {
        expect_type(fun, "closure")
        x <- fun(object)
        expect_s3_class(x, "ggplot")
        x <- fun(object, interestingGroups = "sampleName")
        expect_s3_class(x, "ggplot")
    }
})
