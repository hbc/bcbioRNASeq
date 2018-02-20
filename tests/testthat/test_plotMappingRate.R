context("plotMappingRate")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("plotMappingRate", {
    p <- plotMappingRate(bcb)
    expect_is(p, "ggplot")
})
