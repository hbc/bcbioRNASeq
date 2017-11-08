context("plotTotalReads")

test_that("plotTotalReads", {
    p <- plotTotalReads(bcb)
    expect_is(p, "ggplot")
})
