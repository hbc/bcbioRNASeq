context("subset")

bcb <- examples[["bcb"]]

test_that("sample selection", {
    subset <- bcb[, 1:2]
    expect_equal(
        dim(subset),
        c(505, 2)
    )
    expect_equal(
        colnames(subset),
        c("group1_1", "group1_2")
    )
    # Check that subsetting by name works
    expect_equal(
        bcb[, c("group1_1", "group1_2")],
        subset
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

test_that("gene selection", {
    subset <- bcb[1:3, ]
    expect_equal(
        rownames(subset),
        c("ENSMUSG00000002459",
          "ENSMUSG00000004768",
          "ENSMUSG00000005886")
    )
    # Check that subsetting by name works
    expect_equal(
        subset,
        bcb[c("ENSMUSG00000002459",
              "ENSMUSG00000004768",
              "ENSMUSG00000005886"), ]
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
    # Selecting fewer than 3 genes will generate a warning
    expect_warning(
        bcb[1:2, ],
        "Estimated rdf < 1.0; not estimating variance"
    )
})

test_that("invalid ranges", {
    # Improve the error checking in our subset method. This comes from DESeq2.
    expect_error(
        bcb[0, 0],
        "all samples have 0 counts for all genes. check the counting script."
    )
    # Also improve the error message for single gene selection.
    # Currently too cryptic.
    expect_error(
        bcb[1, ],
        "invalid class \"SummarizedExperiment\" object"
    )
    # Single sample selection fails because of the internal DESeqDataSet
    expect_error(
        suppressWarnings(bcb[, 1]),
        "all genes have equal values for all samples."
    )
})
