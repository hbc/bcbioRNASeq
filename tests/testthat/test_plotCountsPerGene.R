context("plotCountsPerGene")

test_that("plotCountsPerGene", {
    p <- plotCountsPerGene(bcb)
    expect_true(is(p, "ggplot"))
})
