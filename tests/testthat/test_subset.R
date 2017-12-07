context("subset")

bcb <- examples[["bcb"]]

test_that("Normal gene and sample selection", {
    subset <- suppressMessages(bcb[1:100, 1:4])
    expect_equal(
        dim(subset),
        c(100, 4)
    )
    expect_equal(
        rownames(subset)[[1]],
        "ENSMUSG00000002459"
    )
    expect_equal(
        colnames(subset),
        c("group1_1",
          "group1_2",
          "group2_1",
          "group2_2")
    )
    # Check that the internal DESeqDataSet gets updated
    expect_equal(
        dim(subset),
        dim(bcbio(subset, "DESeqDataSet"))
    )
    expect_equal(
        colnames(subset),
        colnames(bcbio(subset, "DESeqDataSet"))
    )
})

test_that("Minimal sample selection", {
    subset <- suppressMessages(bcb[, 1:2])
    expect_equal(
        dim(subset),
        c(505, 2)
    )
    # Check that subsetting by name also works
    expect_equal(
        suppressMessages(bcb[, c("group1_1", "group1_2")]),
        subset
    )
})

test_that("Minimal gene selection", {
    subset <- suppressMessages(
        bcb[1:2, , skipNorm = TRUE]
    )
    expect_equal(
        rownames(subset),
        c("ENSMUSG00000002459",
          "ENSMUSG00000004768")
    )
    # Check that the internal DESeqDataSet gets updated
    expect_equal(
        dim(subset),
        dim(bcbio(subset, "DESeqDataSet"))
    )
    expect_equal(
        rownames(subset),
        rownames(bcbio(subset, "DESeqDataSet"))
    )
    # Check that subsetting by name also works
    expect_equal(
        suppressWarnings(suppressMessages(
            bcb[c("ENSMUSG00000002459",
                  "ENSMUSG00000004768"), , skipNorm = TRUE]
        )),
        subset
    )
    # Selecting fewer than 3 genes will generate a warning
    expect_warning(
        suppressMessages(bcb[1:2, ]),
        "Estimated rdf < 1.0; not estimating variance"
    )
})

test_that("Invalid ranges", {
    expect_error(
        bcb[1, ],
        "At least 2 genes are required"
    )
    expect_error(
        bcb[, 1],
        "At least 2 samples are required"
    )
})
