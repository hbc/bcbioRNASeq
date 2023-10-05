test_that("plotFeatureSaturation", {
    x <- plotFeatureSaturation(
        object = object,
        trendline = TRUE,
        label = TRUE
    )
    expect_s3_class(x, "ggplot")
})
