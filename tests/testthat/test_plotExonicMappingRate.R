context("plotExonicMappingRate")

test_that("plotExonicMappingRate", {
    p <- plotExonicMappingRate(bcb)
    expect_is(p, "ggplot")
})
