context("plotCountsPerGene")

test_that("plotCountsPerGene", {
    p <- plotCountsPerGene(bcb)
    expect_is(p, "ggplot")
})
