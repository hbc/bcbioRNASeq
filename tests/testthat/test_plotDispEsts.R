context("plotDispEsts")

test_that("plotDispEsts", {
    p <- plotDispEsts(bcb)
    expect_true(is(p, "list"))
    expect_equal(
        names(p),
        c("rect", "text")
    )
})
