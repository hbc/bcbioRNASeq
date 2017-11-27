context("plotTotalReads")

bcb <- examples[["bcb"]]

test_that("plotTotalReads", {
    p <- plotTotalReads(bcb)
    expect_is(p, "ggplot")
})
