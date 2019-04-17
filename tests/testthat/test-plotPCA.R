context("plotPCA")

test_that("label", {
    x <- plotPCA(object, label = FALSE)
    expect_s3_class(x, "ggplot")
    x <- plotPCA(object, label = TRUE)
    expect_s3_class(x, "ggplot")
})

test_that("return", {
    x <- plotPCA(object, return = "DataFrame")
    expect_s4_class(x, "DataFrame")
})
