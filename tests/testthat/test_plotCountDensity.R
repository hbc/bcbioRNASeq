context("plotCountDensity")

test_that("plotCountDensity", {
    p <- plotCountDensity(bcb)
    expect_true(is(p, "ggplot"))
})
