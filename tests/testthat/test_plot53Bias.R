context("plot53Bias")

test_that("plot53Bias", {
    p <- plot53Bias(bcb)
    expect_is(p, "ggplot")
})
