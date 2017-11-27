context("plotPCA")

bcb <- examples[["bcb"]]

test_that("plotPCA", {
    p <- plotPCA(bcb)
    expect_is(p, "ggplot")
})
