context("plotDispEsts")

## NOTE This can return "figure margins too large" error.

test_that("plotDispEsts", {
    x <- plotDispEsts(object)
    expect_type(x, "list")
    expect_identical(
        object = names(x),
        expected = c("rect", "text")
    )
})
