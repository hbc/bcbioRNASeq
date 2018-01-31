context("internal")

test_that("tx2gene", {
    uploadDir <- system.file(
        file.path("extdata", "bcbio"),
        package = "bcbioRNASeq")
    file <- file.path(uploadDir, "2017-05-23_rnaseq")
    df <- .tx2gene(file, organism = "hg19")
    expect_equal(ncol(df), 2L)
})
