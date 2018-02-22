context("subset")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("Normal gene and sample selection", {
    subset <- suppressMessages(bcb[1L:100L, 1L:4L])
    expect_identical(
        dim(subset),
        c(100L, 4L)
    )
    expect_identical(
        rownames(subset)[[1L]],
        "ENSMUSG00000002459"
    )
    expect_identical(
        colnames(subset),
        c("group1_1",
          "group1_2",
          "group2_1",
          "group2_2")
    )
    # Check that the internal DESeqDataSet gets updated
    expect_identical(
        dim(subset),
        dim(bcbio(subset, "DESeqDataSet"))
    )
    expect_identical(
        colnames(subset),
        colnames(bcbio(subset, "DESeqDataSet"))
    )
})

test_that("Minimal sample selection", {
    subset <- suppressMessages(bcb[, 1L:2L])
    expect_identical(
        dim(subset),
        c(505L, 2L)
    )
    # Check that subsetting by name also works
    expect_identical(
        suppressMessages(bcb[, c("group1_1", "group1_2")]),
        subset
    )
})

test_that("Minimal gene selection", {
    subset <- suppressWarnings(suppressMessages(
        bcb[1L:2L, , transform = FALSE]
    ))
    expect_identical(
        rownames(subset),
        c("ENSMUSG00000002459",
          "ENSMUSG00000004768")
    )
    # Check that the internal DESeqDataSet gets updated
    expect_identical(
        dim(subset),
        dim(bcbio(subset, "DESeqDataSet"))
    )
    expect_identical(
        rownames(subset),
        rownames(bcbio(subset, "DESeqDataSet"))
    )
    # Check that subsetting by name also works
    expect_identical(
        suppressWarnings(suppressMessages(
            bcb[c("ENSMUSG00000002459",
                  "ENSMUSG00000004768"), , transform = FALSE]
        )),
        subset
    )
    # Selecting fewer than 3 genes will generate a warning
    expect_warning(
        suppressMessages(bcb[1L:2L, ]),
        "Estimated rdf < 1.0; not estimating variance"
    )
})

test_that("Invalid ranges", {
    expect_error(
        bcb[1L, ],
        "is_greater_than : length\\(i\\) are not all greater than 1L."
    )
    expect_error(
        bcb[, 1L],
        "is_greater_than : length\\(j\\) are not all greater than 1L."
    )
})
