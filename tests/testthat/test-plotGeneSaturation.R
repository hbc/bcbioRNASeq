context("plotGeneSaturation")

test_that("plotGeneSaturation", {
    x <- plotGeneSaturation(
        object = object,
        trendline = TRUE,
        label = TRUE
    )
    expect_s3_class(x, "ggplot")
})
