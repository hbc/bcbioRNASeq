context("plotIntronicMappingRate")

test_that("plotIntronicMappingRate", {
    p <- plotIntronicMappingRate(bcb_small)
    expect_is(p, "ggplot")
})
