context("plotRRNAMappingRate")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("plotRRNAMappingRate", {
    p <- plotRRNAMappingRate(bcb)
    expect_is(p, "ggplot")
})
