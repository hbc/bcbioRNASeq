context("plotFeaturesDetected")

test_that("bcbioRNASeq", {
    x <- plotFeaturesDetected(bcb)
    expect_s3_class(x, "ggplot")
})



context("plotGenesDetected")

test_that("bcbioRNASeq", {
    x <- plotGenesDetected(bcb)
    expect_s3_class(x, "ggplot")
})
