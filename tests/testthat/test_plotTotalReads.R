context("plotTotalReads")

test_that("plotTotalReads", {
    p <- plotTotalReads(bcb_small)
    expect_is(p, "ggplot")
})
