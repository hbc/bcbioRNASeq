context("plotGeneSaturation")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("plotGeneSaturation", {
    p <- plotGeneSaturation(bcb)
    expect_is(p, "ggplot")
})
