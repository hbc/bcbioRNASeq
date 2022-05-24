test_that("Clean mode", {
    x <- sampleData(object, clean = TRUE)
    ## Require that all clean columns are factor.
    invisible(lapply(x, function(x) {
        expect_is(x, "factor")
    }))
})

test_that("Verbose mode", {
    x <- sampleData(object, clean = FALSE)
    ## Return `interestingGroups` factor column by default.
    expect_is(x[["interestingGroups"]], "factor")
    ## Otherwise it should be identical to `colData`.
    x[["interestingGroups"]] <- NULL
    expect_identical(
        object = as.data.frame(x),
        expected = as.data.frame(camelCase(colData(object), strict = TRUE))
    )
})

test_that("sampleData<-", {
    sampleData(object)[["testthat"]] <- factor("XXX")
    expect_identical(
        object = levels(sampleData(object)[["testthat"]]),
        expected = "XXX"
    )
})
