context("Internal Functions")

test_that("tx2gene", {
    df <- .tx2gene(
        projectDir = file.path(uploadDir, "2017-05-23_rnaseq")
    )
    expect_identical(ncol(df), 2L)
})
