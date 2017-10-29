context("plotMappingRate")

test_that("plotMappingRate", {
    p <- plotMappingRate(bcb)
    expect_is(p, "ggplot")
})
