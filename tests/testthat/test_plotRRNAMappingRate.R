context("plotRRNAMappingRate")

test_that("plotRRNAMappingRate", {
    p <- plotRRNAMappingRate(bcb)
    expect_true(is(p, "ggplot"))
})
