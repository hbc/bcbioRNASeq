context("plotCountsPerGene")

bcb <- examples[["bcb"]]

test_that("plotCountsPerGene", {
    p <- plotCountsPerGene(bcb)
    expect_is(p, "ggplot")
})
