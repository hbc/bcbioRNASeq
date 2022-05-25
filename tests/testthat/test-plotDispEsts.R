## NOTE This can return "figure margins too large" error.

test_that("plotDispEsts", {
    x <- plotDispEsts(object)
    expect_type(x, "list")
    expect_named(
        object = x,
        expected = c("rect", "text")
    )
})
