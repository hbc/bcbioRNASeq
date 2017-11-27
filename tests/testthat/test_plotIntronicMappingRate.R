context("plotIntronicMappingRate")

bcb <- examples[["bcb"]]

test_that("plotIntronicMappingRate", {
    p <- plotIntronicMappingRate(bcb)
    expect_is(p, "ggplot")
})
