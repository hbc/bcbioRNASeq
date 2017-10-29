context("plotMappingRate")

test_that("plotMappingRate", {
    p <- plotMappingRate(bcb)
    expect_true(is(p, "ggplot"))
})
