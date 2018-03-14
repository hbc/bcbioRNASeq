context("counts")

test_that("normalized", {
    normalized <- list(FALSE, TRUE, "tpm", "tmm", "rlog", "vst")

    # Check that all are matrices
    expect_true(all(vapply(
        X = normalized,
        FUN = function(arg) {
            is.matrix(counts(bcb_small, normalized = arg))
        },
        FUN.VALUE = logical(1L)
    )))

    # FALSE
    expect_identical(
        counts(bcb_small, normalized = FALSE),
        assays(bcb_small)[["raw"]]
    )
    expect_identical(
        counts(bcb_small, normalized = FALSE),
        assay(bcb_small)
    )

    # TRUE = calculated on the fly with DESeq2

    # tpm
    expect_identical(
        counts(bcb_small, normalized = "tpm"),
        assays(bcb_small)[["tpm"]]
    )
    expect_identical(
        counts(bcb_small, normalized = "tpm"),
        tpm(bcb_small)
    )

    # tmm = calculated on the fly
    expect_identical(
        counts(bcb_small, normalized = "tmm"),
        tmm(bcb_small)
    )

    # rlog
    expect_identical(
        counts(bcb_small, normalized = "rlog"),
        assays(bcb_small)[["rlog"]]
    )

    # vst
    expect_identical(
        counts(bcb_small, normalized = "vst"),
        assays(bcb_small)[["vst"]]
    )
})

test_that("transformationLimit : skipped DESeq transforms", {
    skip <- bcb_small
    # Using `assays<-` will coerce bcbioRNASeq to SummarizedExperiment
    slot(skip, "assays")[["rlog"]] <- NULL
    expect_warning(
        counts(skip, normalized = "rlog"),
        paste(
            "rlog not present in assays.",
            "Calculating log2 TMM counts instead."
        )
    )
    counts <- suppressWarnings(counts(skip, normalized = "rlog"))
    expect_is(counts, "matrix")
})
