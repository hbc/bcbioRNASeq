context("plotRRNAMappingRate")

test_that("plotRRNAMappingRate", {
    p <- plotRRNAMappingRate(bcb)
    expect_is(p, "ggplot")
})
