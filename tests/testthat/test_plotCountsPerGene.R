context("plotCountsPerGene")

test_that("plotCountsPerGene", {
    p <- plotCountsPerGene(bcb_small)
    expect_is(p, "ggplot")
})
