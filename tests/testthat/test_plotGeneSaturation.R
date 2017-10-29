context("plotGeneSaturation")

test_that("plotGeneSaturation", {
    p <- plotGeneSaturation(bcb)
    expect_true(is(p, "ggplot"))
})
