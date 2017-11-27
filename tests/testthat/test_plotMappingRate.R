context("plotMappingRate")

bcb <- examples[["bcb"]]

test_that("plotMappingRate", {
    p <- plotMappingRate(bcb)
    expect_is(p, "ggplot")
})
