context("plotPCA")

test_that("label", {
    x <- plotPCA(object, label = FALSE)
    expect_s3_class(x, "ggplot")
    x <- plotPCA(object, label = TRUE)
    expect_s3_class(x, "ggplot")
})

test_that("Fast mode", {
    expect_error(
        object = plotPCA(bcb_fast),
        expected = "fast mode"
    )
})
