context("relativeLogExpression")

test_that("SummarizedExperiment", {
    x <- relativeLogExpression(object)
    expect_is(x, "matrix")
    expect_identical(
        object = round(unname(x[1L, , drop = TRUE]), digits = 3L),
        expected = c(1.388, 1.140, 1.266, 0.822, 0.916, 1.313)
    )
})
