context("plot53Bias")

test_that("plot53Bias", {
    p <- plot53Bias(bcb_small)
    expect_is(p, "ggplot")
})
