context("sampleData")

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
    ## This check fails due to "DFrame"/"DataFrame" class mismatch on
    ## Bioconductor 3.10, unless we coerce to DataFrame here.
    expect_identical(
        object = as(x, "DataFrame"),
        expected = colData(object)
    )
})

test_that("sampleData<-", {
    sampleData(object)[["testthat"]] <- factor("XXX")
    expect_identical(
        object = levels(sampleData(object)[["testthat"]]),
        expected = "XXX"
    )
})
