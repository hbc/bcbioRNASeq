context("flatFiles")

bcb <- examples[["bcb"]]

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
