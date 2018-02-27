context("flatFiles")

test_that("flatFiles", {
    flat <- flatFiles(bcb)
    expect_is(flat, "list")
    expect_identical(
        names(flat),
        c(
            "assays",
            "rowData",
            "colData",
            "metadata",
            "bcbio"
        )
    )
})
