context("sampleData")

test_that("Clean mode", {
    object <- bcb
    x <- sampleData(object, clean = TRUE)
    # Require that all clean columns are factor.
    invisible(lapply(x, function(x) {
        expect_is(x, "factor")
    }))
})

test_that("Verbose mode", {
    object <- bcb
    x <- sampleData(object, clean = FALSE)

    # Return `interestingGroups` factor column by default.
    expect_is(x[["interestingGroups"]], "factor")

    # Otherwise it should be identical to `colData`.
    x[["interestingGroups"]] <- NULL
    expect_identical(
        object = x,
        expected = colData(object)
    )
})

test_that("sampleData<-", {
    object <- bcb
    sampleData(object)[["testthat"]] <- factor("XXX")
    expect_identical(
        object = levels(sampleData(object)[["testthat"]]),
        expected = "XXX"
    )
})
