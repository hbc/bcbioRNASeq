context("Internal Functions")

test_that("tx2gene", {
    file <- file.path(uploadDir, "2017-05-23_rnaseq")
    df <- .tx2gene(file, organism = "hg19")
    expect_identical(ncol(df), 2L)
})
