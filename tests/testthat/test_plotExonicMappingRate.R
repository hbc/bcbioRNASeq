context("plotExonicMappingRate")

test_that("plotExonicMappingRate", {
    p <- plotExonicMappingRate(bcb)
    expect_true(is(p, "ggplot"))
})
