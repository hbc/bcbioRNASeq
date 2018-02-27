context("plotDispEsts")

test_that("plotDispEsts", {
    p <- plotDispEsts(bcb)
    expect_is(p, "list")
    expect_identical(
        names(p),
        c("rect", "text")
    )
})
