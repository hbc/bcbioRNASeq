context("plotRRNAMappingRate")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("plotRRNAMappingRate", {
    p <- plotRRNAMappingRate(bcb)
    expect_is(p, "ggplot")
})
