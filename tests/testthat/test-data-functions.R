context("Data Functions")

bcb <- bcbioRNASeq(
    uploadDir = system.file("extdata/bcbio", package = "bcbioRNASeq"),
    rlog = TRUE,
    vst = TRUE
)
dds <- as(bcb, "DESeqDataSet")



# aggregateReplicates ==========================================================
test_that("aggregateReplicates", {
    # Assign groupings into `aggregate` column of `colData()`
    aggregate <- as.factor(sub("^([a-z0-9]+)_.*", "\\1", colnames(bcb)))
    names(aggregate) <- colnames(bcb)
    bcb[["aggregate"]] <- aggregate
    x <- aggregateReplicates(bcb)
    expect_identical(colnames(x), c("group1", "group2"))
    expect_identical(sum(counts(x)), sum(counts(bcb)))
    expect_equal(rowSums(counts(x)), rowSums(counts(bcb)))
})



# counts =======================================================================
test_that("counts : normalized argument", {
    normalized <- list(FALSE, TRUE, "tpm", "tmm", "vst")

    # Check that all are matrices
    expect_true(all(vapply(
        X = normalized,
        FUN = function(arg) {
            is.matrix(counts(bcb, normalized = arg))
        },
        FUN.VALUE = logical(1L)
    )))

    # FALSE
    expect_identical(
        counts(bcb, normalized = FALSE),
        assays(bcb)[["counts"]]
    )
    expect_identical(
        counts(bcb, normalized = FALSE),
        assay(bcb)
    )

    # TRUE
    expect_identical(
        counts(bcb, normalized = TRUE),
        assays(bcb)[["normalized"]]
    )

    # tpm
    expect_identical(
        counts(bcb, normalized = "tpm"),
        assays(bcb)[["tpm"]]
    )
    expect_identical(
        counts(bcb, normalized = "tpm"),
        tpm(bcb)
    )

    # tmm: calculated on the fly
    expect_identical(
        counts(bcb, normalized = "tmm"),
        tmm(bcb)
    )

    # rle: calculated on the fly
    expect_is(
        counts(bcb, normalized = "rle"),
        "matrix"
    )

    # rlog
    expect_identical(
        counts(bcb, normalized = "rlog"),
        assays(bcb)[["rlog"]]
    )

    # vst
    expect_identical(
        counts(bcb, normalized = "vst"),
        assays(bcb)[["vst"]]
    )
})

test_that("counts : skipped DESeq2 transforms", {
    skip <- bcb
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



# tmm ==========================================================================
test_that("tmm", {
    expect_is(tmm(bcb), "matrix")
    expect_is(tmm(dds), "matrix")
    expect_is(tmm(assay(bcb)), "matrix")
})
