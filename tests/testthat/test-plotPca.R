test_that("label", {
    x <- plotPca(object, label = FALSE)
    expect_s3_class(x, "ggplot")
    x <- plotPca(object, label = TRUE)
    expect_s3_class(x, "ggplot")
})

test_that("Fast mode", {
    expect_error(
        object = plotPca(bcb_fast),
        regexp = "fast mode"
    )
})
