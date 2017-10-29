context("plotTotalReads")

test_that("plotTotalReads", {
    p <- plotTotalReads(bcb)
    expect_true(is(p, "ggplot"))
})
