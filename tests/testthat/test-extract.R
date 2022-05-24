nrow <- 50L
ncol <- 2L

test_that("Normal gene and sample selection", {
    subset <- object[seq_len(nrow), seq_len(ncol)]
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        object = dim(subset),
        expected = c(nrow, ncol)
    )
    expect_identical(
        object = rownames(subset),
        expected = head(rownames(object), n = nrow)
    )
    expect_identical(
        object = colnames(subset),
        expected = head(colnames(object), n = ncol)
    )
})

test_that("Require at least 50 genes, 2 samples", {
    expect_error(object[seq_len(49L), ])
    expect_error(object[, seq_len(1L)])
})

test_that("Check for unmodified return when using empty brackets", {
    expect_identical(
        object = object[, ],
        expected = object
    )
})

test_that("Calculate DESeq2 transforms by default", {
    ## Transform enabled by default (if calculated).
    expect_identical(
        object = sort(assayNames(object[seq_len(nrow), ])),
        expected = c(
            "aligned",
            "avgTxLength",
            "counts",
            "fpkm",
            "normalized",
            "tpm",
            "vst"
        )
    )
})

test_that("Allow the user to skip transforms, using 'recalculate'", {
    expect_identical(
        object = sort(assayNames(object[seq_len(nrow), , recalculate = FALSE])),
        expected = c(
            "aligned",
            "avgTxLength",
            "counts",
            "tpm"
        )
    )
})
