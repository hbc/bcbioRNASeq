context("plotIntronicMappingRate")

test_that("plotIntronicMappingRate", {
    p <- plotIntronicMappingRate(bcb)
    expect_true(is(p, "ggplot"))
})
