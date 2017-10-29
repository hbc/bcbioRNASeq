context("plotPCA")

test_that("plotPCA", {
    p <- plotPCA(bcb)
    expect_is(p, "ggplot")
})
