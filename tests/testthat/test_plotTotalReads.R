context("plotTotalReads")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("plotTotalReads", {
    p <- plotTotalReads(bcb)
    expect_is(p, "ggplot")
})
