skip_if_not_installed("vsn")

test_that("plotMeanSd", {
    x <- plotMeanSd(object)
    expect_s3_class(x, "ggplot")
})
