context("plotMappedReads")

bcb <- examples[["bcb"]]

test_that("plotMappedReads", {
    p <- plotMappedReads(bcb)
    expect_is(p, "ggplot")
})
