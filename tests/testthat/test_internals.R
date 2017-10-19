context("Internals")

test_that("tx2gene", {
    dataDir <- system.file(file.path("extdata", "bcbio"),
                             package = "bcbioRNASeq") %>%
        normalizePath()
    fn <- file.path(dataDir, "2017-05-23_rnaseq")
    df <- .tx2gene(fn, "hg19")
    expect_equal(ncol(df), 2)
})
