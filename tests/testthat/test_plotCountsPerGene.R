context("plotCountsPerGene")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("plotCountsPerGene", {
    p <- plotCountsPerGene(bcb)
    expect_is(p, "ggplot")
})
