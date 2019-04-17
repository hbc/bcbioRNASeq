context("plotCountsPerBiotype")

test_that("plotCountsPerBiotype", {
    x <- plotCountsPerBiotype(object)
    expect_s3_class(x, "ggplot")
})



context("plotCountsPerBroadClass")

test_that("plotCountsPerBroadClass", {
    x <- plotCountsPerBroadClass(object)
    expect_s3_class(x, "ggplot")
})
