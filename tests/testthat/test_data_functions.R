context("Data Functions")



# aggregateReplicates ==========================================================
test_that("aggregateReplicates", {
    bcb <- bcb_small
    # Assign groupings into `aggregate` column of `colData()`
    aggregate <- as.factor(sub("^([a-z]+)_.*", "\\1", colnames(bcb)))
    names(aggregate) <- colnames(bcb)
    bcb[["aggregate"]] <- aggregate
    x <- aggregateReplicates(bcb)
    expect_identical(dim(x), c(500L, 2L))
    expect_identical(sum(counts(x)), sum(counts(bcb)))
    expect_equal(rowSums(counts(x)), rowSums(counts(bcb)))
})



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
        assays(bcb_small)[["counts"]]
    )
    expect_identical(
        counts(bcb_small, normalized = FALSE),
        assay(bcb_small)
    )

    # TRUE
    expect_identical(
        counts(bcb_small, normalized = TRUE),
        assays(bcb_small)[["normalized"]]
    )

    # tpm
    expect_identical(
        counts(bcb_small, normalized = "tpm"),
        assays(bcb_small)[["tpm"]]
    )
    expect_identical(
        counts(bcb_small, normalized = "tpm"),
        tpm(bcb_small)
    )

    # tmm: calculated on the fly
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



# sampleData ===================================================================
# Check output of `return` parameter
return <- methodFormals(
    f = "sampleData",
    signature = "bcbioRNASeq"
) %>%
    .[["return"]] %>%
    as.character() %>%
    .[-1L]

test_that("sampleData: Verbose mode", {
    list <- lapply(return, function(x) {
        sampleData(
            bcb_small,
            clean = FALSE,
            return = x
        )
    })

    # Check returns
    expect_identical(
        lapply(list, class),
        list(
            structure("DataFrame", package = "S4Vectors"),
            "data.frame",
            "knitr_kable"
        )
    )

    # Check dimnames
    expected <- list(
        colnames(bcb_small),
        c(colnames(colData(bcb_small)), "interestingGroups")
    )
    expect_identical(
        lapply(list, dimnames),
        list(
            expected,
            expected,
            NULL
        )
    )
})

test_that("sampleData : Clean mode", {
    list <- lapply(return, function(x) {
        sampleData(
            bcb_small,
            clean = TRUE,
            return = x
        )
    })

    # Check returns
    expect_identical(
        lapply(list, class),
        list(
            structure("DataFrame", package = "S4Vectors"),
            "data.frame",
            "knitr_kable"
        )
    )

    # Check dimnames
    expected <- list(
        colnames(bcb_small),
        c(
            "sampleName",
            "day",
            "replicate",
            "strain",
            "tissue",
            "treatment"
        )
    )
    expect_identical(
        lapply(list, dimnames),
        list(
            expected,
            expected,
            NULL
        )
    )

    # Ensure all columns are factor
    invisible(lapply(list[[1L]], function(x) {
        expect_is(x, "factor")
    }))

    # Interesting groups
    x <- sampleData(bcb_small, clean = FALSE, interestingGroups = NULL)
    expect_identical(
        x[["interestingGruops"]],
        NULL
    )
    x <- sampleData(bcb_small, clean = FALSE, interestingGroups = "day")
    expect_identical(
        levels(x[["interestingGroups"]]),
        c("0", "7")
    )
})



# selectSamples ================================================================
test_that("selectSamples : bcbioRNASeq", {
    x <- selectSamples(bcb_small, treatment = "folic_acid")
    expect_identical(dim(x), c(500L, 3L))
    expect_identical(
        names(assays(x)),
        c("counts", "tpm", "length", "normalized", "rlog", "vst")
    )
})

test_that("selectSamples : DESeqDataSet", {
    x <- selectSamples(dds_small, treatment = "folic_acid")
    expect_identical(dim(x), c(500L, 3L))
})
