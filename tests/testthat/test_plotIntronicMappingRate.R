context("plotIntronicMappingRate")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("plotIntronicMappingRate", {
    p <- plotIntronicMappingRate(bcb)
    expect_is(p, "ggplot")
})
