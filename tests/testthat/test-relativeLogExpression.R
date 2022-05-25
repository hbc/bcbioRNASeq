test_that("SummarizedExperiment", {
    x <- relativeLogExpression(object)
    expect_type(x, "double")
    expect_identical(
        object = round(unname(x[1L, , drop = TRUE]), digits = 2L),
        expected = c(29079.47, 39937.62, 28051.07, 33200.52, 33200.52, 35928.23)
    )
})
