context("plot53Bias")

bcb <- examples[["bcb"]]

test_that("plot53Bias", {
    p <- plot53Bias(bcb)
    expect_is(p, "ggplot")
})
