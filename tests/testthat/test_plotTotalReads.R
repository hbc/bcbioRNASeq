context("plotTotalReads")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("plotTotalReads", {
    p <- plotTotalReads(bcb)
    expect_is(p, "ggplot")
})
