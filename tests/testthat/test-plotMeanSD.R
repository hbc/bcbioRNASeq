context("plotMeanSD")

test_that("plotMeanSD", {
    x <- plotMeanSD(object)
    expect_s3_class(x, "ggplot")
})
