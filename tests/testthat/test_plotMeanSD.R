context("plotMeanSD")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("plotGenderMarkers", {
    p <- plotMeanSD(bcb)
    expect_is(p, "ggplot")
})
