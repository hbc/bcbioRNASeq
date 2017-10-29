context("plotCountDensity")

test_that("plotCountDensity", {
    p <- plotCountDensity(bcb)
    expect_is(p, "ggplot")
})
