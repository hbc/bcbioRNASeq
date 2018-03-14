context("plotGeneSaturation")

test_that("plotGeneSaturation", {
    p <- plotGeneSaturation(bcb_small)
    expect_is(p, "ggplot")
})
