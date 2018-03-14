context("plotCountDensity")

test_that("plotCountDensity", {
    p <- plotCountDensity(bcb_small)
    expect_is(p, "ggplot")
})
