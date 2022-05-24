skip_if_not_installed("vsn")

test_that("plotMeanSD", {
    x <- plotMeanSD(object)
    expect_s3_class(x, "ggplot")
})
