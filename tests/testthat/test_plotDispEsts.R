context("plotDispEsts")

test_that("plotDispEsts", {
    p <- plotDispEsts(bcb_small)
    expect_is(p, "list")
    expect_identical(
        names(p),
        c("rect", "text")
    )
})
