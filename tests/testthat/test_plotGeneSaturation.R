context("plotGeneSaturation")

test_that("plotGeneSaturation", {
    p <- plotGeneSaturation(bcb)
    expect_is(p, "ggplot")
})
