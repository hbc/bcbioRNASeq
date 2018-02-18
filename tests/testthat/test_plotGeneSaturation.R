context("plotGeneSaturation")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("plotGeneSaturation", {
    p <- plotGeneSaturation(bcb)
    expect_is(p, "ggplot")
})
