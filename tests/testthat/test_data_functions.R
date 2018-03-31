context("Data Functions")



# counts =======================================================================
test_that("counts : normalized argument", {
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

test_that("counts : apply transformationLimit", {
    skip <- bcb_small
    # Using `assays<-` will coerce bcbioRNASeq to SummarizedExperiment
    slot(skip, "assays")[["rlog"]] <- NULL
    slot(skip, "assays")[["vst"]] <- NULL
    expect_warning(
        counts(skip, normalized = "rlog"),
        paste(
            "rlog not present in assays.",
            "Calculating log2 TMM counts instead."
        )
    )
    expect_warning(
        counts(skip, normalized = "vst"),
        paste(
            "vst not present in assays.",
            "Calculating log2 TMM counts instead."
        )
    )
    counts <- suppressWarnings(counts(skip, normalized = "rlog"))
    expect_is(counts, "matrix")
})



# interestingGroups ============================================================
test_that("interestingGroups", {
    expect_identical(
        interestingGroups(bcb_small),
        "treatment"
    )
})

test_that("interestingGroups<-", {
    expect_silent(
        interestingGroups(bcb_small) <- "sampleName"
    )
    # Interesting group must be in the metadata
    expect_error(
        interestingGroups(bcb_small) <- "XXX",
        paste(
            "is_subset :",
            "The element 'XXX' in interestingGroups is not in",
            "colnames\\(x\\)."
        )
    )
})



# selectSamples ================================================================
test_that("selectSamples : bcbioRNASeq", {
    x <- selectSamples(
        object = bcb_small,
        treatment = "folic_acid",
        transform = TRUE
    )
    expect_identical(dim(x), c(500L, 3L))
    expect_identical(
        names(assays(x)),
        c("raw", "tpm", "length", "rlog", "vst")
    )
    x <- selectSamples(
        object = bcb_small,
        treatment = "folic_acid",
        transform = FALSE
    )
    expect_identical(dim(x), c(500L, 3L))
    expect_identical(
        names(assays(x)),
        c("raw", "tpm", "length")
    )
})

test_that("selectSamples : DESeqDataSet", {
    x <- selectSamples(dds_small, treatment = "folic_acid")
    expect_identical(dim(x), c(500L, 3L))
})
