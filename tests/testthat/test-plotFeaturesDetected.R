test_that("bcbioRNASeq", {
    x <- plotFeaturesDetected(bcb)
    expect_s3_class(x, "ggplot")
})
