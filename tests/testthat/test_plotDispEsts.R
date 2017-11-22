context("plotDispEsts")

bcb <- examples[["bcb"]]

test_that("plotDispEsts", {
    p <- plotDispEsts(bcb)
    expect_is(p, "list")
    expect_equal(
        names(p),
        c("rect", "text")
    )
})
