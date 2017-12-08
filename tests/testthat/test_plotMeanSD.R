context("plotMeanSD")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("plotGenderMarkers", {
    p <- plotMeanSD(bcb)
    expect_is(p, "ggplot")
})
