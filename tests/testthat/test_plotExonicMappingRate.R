context("plotExonicMappingRate")

test_that("plotExonicMappingRate", {
    p <- plotExonicMappingRate(bcb_small)
    expect_is(p, "ggplot")
})
