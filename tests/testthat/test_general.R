context("General")

data(bcb, envir = environment())
assay <- SummarizedExperiment::assay



# counts =======================================================================
with_parameters_test_that(
    "counts : Slotted assays", {
        object <- bcb
        # Check that all are matrices.
        expect_is(
            object = counts(object, normalized = normalized),
            class = "matrix"
        )
        # Check that we're matching the expected assay matrix.
        expect_identical(
            object = counts(object, normalized = normalized),
            expected = assays(object)[[assay]]
        )
    },
    normalized = list(
        FALSE,
        TRUE,
        "tpm",
        "vst"
    ),
    assay = list(
        "counts",
        "normalized",
        "tpm",
        "vst"
    )
)

with_parameters_test_that(
    "counts : On the fly assays", {
        object <- bcb
        expect_is(
            object = counts(object, normalized = normalized),
            class = "matrix"
        )
        expect_null(assays(object)[[normalized]])
    },
    normalized = c("tmm", "rle")
)

with_parameters_test_that(
    "counts : Skipped DESeq2 transforms", {
        object <- bcb
        assays(object)[[normalized]] <- NULL
        expect_error(
            object = counts(object, normalized = normalized),
            regexp = "assayNames"
        )
    },
    normalized = c("rlog", "vst")
)



# sampleData ===================================================================
test_that("sampleData : Verbose mode (default)", {
    object <- bcb
    x <- sampleData(object, clean = FALSE)

    # Return `interestingGroups` factor column by default.
    expect_is(x[["interestingGroups"]], "factor")

    # Otherwise it should be identical to `colData`.
    x[["interestingGroups"]] <- NULL
    expect_identical(
        object = x,
        expected = colData(object)
    )
})

test_that("sampleData : Clean mode", {
    object <- bcb
    x <- sampleData(object, clean = TRUE)
    # Require that all clean columns are factor.
    invisible(lapply(x, function(x) {
        expect_is(x, "factor")
    }))
})

test_that("sampleData<-", {
    object <- bcb
    sampleData(object)[["testthat"]] <- factor("XXX")
    expect_identical(
        object = levels(sampleData(object)[["testthat"]]),
        expected = "XXX"
    )
})



# tmm ==========================================================================
with_parameters_test_that(
    "tmm", {
        expect_is(
            object = tmm(object),
            class = "matrix"
        )
    },
    object = list(
        bcbioRNASeq = bcb,
        matrix = assay(bcb)
    )
)
