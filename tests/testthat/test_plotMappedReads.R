context("plotMappedReads")

test_that("plotMappedReads", {
    p <- plotMappedReads(bcb_small)
    expect_is(p, "ggplot")
})
