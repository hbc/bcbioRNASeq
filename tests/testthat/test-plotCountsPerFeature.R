context("plotCountsPerFeature")

test_that("bcbioRNASeq", {
    x <- plotCountsPerFeature(object, geom = "boxplot")
    expect_s3_class(x, "ggplot")
    x <- plotCountsPerFeature(
        object = object,
        geom = "density",
        normalized = "vst",
        interestingGroups = "sampleName"
    )
    expect_s3_class(x, "ggplot")
})
