context("plotExonicMappingRate")

bcb <- examples[["bcb"]]

test_that("plotExonicMappingRate", {
    p <- plotExonicMappingRate(bcb)
    expect_is(p, "ggplot")
})
