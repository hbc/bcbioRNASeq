context("plotCountDensity")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("plotCountDensity", {
    p <- plotCountDensity(bcb)
    expect_is(p, "ggplot")
})
