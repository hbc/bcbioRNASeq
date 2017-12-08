context("flatFiles")

load(system.file(
    file.path("inst", "extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("flatFiles", {
    flat <- flatFiles(bcb)
    expect_is(flat, "list")
    expect_equal(
        names(flat),
        c("assays",
          "rowData",
          "colData",
          "metadata",
          "bcbio")
    )
})
