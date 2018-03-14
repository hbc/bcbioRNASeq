context("plotMappingRate")

test_that("plotMappingRate", {
    p <- plotMappingRate(bcb_small)
    expect_is(p, "ggplot")
})
