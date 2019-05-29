context("counts")

with_parameters_test_that(
    "Slotted assays", {
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
    "On the fly assay calculations", {
        expect_is(
            object = counts(object, normalized = normalized),
            class = "matrix"
        )
        expect_null(assays(object)[[normalized]])
    },
    normalized = c("tmm", "rle")
)

with_parameters_test_that(
    "Skipped DESeq2 transforms", {
        assays(object)[[normalized]] <- NULL
        expect_error(
            object = counts(object, normalized = normalized),
            regexp = normalized
        )
    },
    normalized = c("rlog", "vst")
)
