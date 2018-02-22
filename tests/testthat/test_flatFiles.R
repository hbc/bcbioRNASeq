context("flatFiles")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("flatFiles", {
    flat <- flatFiles(bcb)
    expect_is(flat, "list")
    expect_identical(
        names(flat),
        c("assays",
          "rowData",
          "colData",
          "metadata",
          "bcbio")
    )
})
