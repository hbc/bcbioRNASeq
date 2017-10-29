context("plotPCA")

test_that("plotPCA", {
    p <- plotPCA(bcb)
    expect_true(is(p, "ggplot"))
})
