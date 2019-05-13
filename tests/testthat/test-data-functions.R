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



# interestingGroups ============================================================
test_that("interestingGroups : NULL handling", {
    expect_is(
        metrics(bcb, interestingGroups = NULL),
        "data.frame"
    )
    expect_is(
        plotTotalReads(bcb, interestingGroups = NULL),
        "ggplot"
    )
    expect_is(
        sampleData(bcb, interestingGroups = NULL),
        "DataFrame"
    )

    x <- bcb
    expect_error(interestingGroups(x) <- NULL)
    metadata(x)[["interestingGroups"]] <- NULL
    expect_error(validObject(x))
})



# sampleData ===================================================================
test_that("sampleData : Verbose mode (default)", {
    # Match colData when `interestingGroups = NULL`
    expect_identical(
        sampleData(bcb, clean = FALSE, interestingGroups = NULL),
        colData(bcb)
    )

    # Return `interestingGroups` factor column by default
    x <- sampleData(bcb, clean = FALSE)
    expect_is(x[["interestingGroups"]], "factor")

    # Interesting groups
    x <- sampleData(bcb, clean = FALSE, interestingGroups = NULL)
    expect_identical(
        x[["interestingGruops"]],
        NULL
    )
    x <- sampleData(bcb, clean = FALSE, interestingGroups = "group")
    expect_identical(
        levels(x[["interestingGroups"]]),
        c("ctrl", "ko")
    )
})

test_that("sampleData : Clean mode", {
    # Return useful factor columns
    x <- sampleData(bcb, clean = TRUE)
    expect_identical(
        colnames(x),
        c("sampleName", "group")
    )
    # Ensure all columns are factor
    invisible(lapply(x, function(x) {
        expect_is(x, "factor")
    }))
})

test_that("sampleData assignment", {
    x <- bcb
    sampleData(x)[["testthat"]] <- factor("XXX")
    expect_identical(levels(sampleData(x)[["testthat"]]), "XXX")
})



# selectSamples ================================================================
test_that("selectSamples : bcbioRNASeq", {
    x <- selectSamples(bcb, group = "ctrl")
    expect_identical(colnames(x), c("group1_1", "group1_2"))
    expect_identical(
        assayNames(x),
        c("counts", "tpm", "length", "normalized", "vst", "rlog")
    )
})

test_that("selectSamples : DESeqDataSet", {
    x <- selectSamples(dds, group = "ctrl")
    expect_identical(colnames(x), c("group1_1", "group1_2"))
})



# tmm ==========================================================================
test_that("tmm", {
    expect_is(tmm(bcb), "matrix")
    expect_is(tmm(dds), "matrix")
    expect_is(tmm(assay(bcb)), "matrix")
})
