context("relativeLogExpression")

test_that("SummarizedExperiment", {
    x <- relativeLogExpression(object)
    expect_is(x, "matrix")
    expect_identical(
        object = round(unname(x[1L, , drop = TRUE]), digits = 3L),
        expected = c(103.579, 89.313, 100.178, 62.151, 71.880, 108.420)
    )
})
