test_that("Slotted assays", {
    Map(
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
        ),
        f = function(normalized, assay) {
            ## Check that all are matrices.
            expect_is(
                object = counts(object, normalized = normalized),
                class = "matrix"
            )
            ## Check that we're matching the expected assay matrix.
            expect_identical(
                object = counts(object, normalized = normalized),
                expected = assays(object)[[assay]]
            )
        }
    )
})

test_that("On the fly assay calculations", {
    for (normalized in c("tmm", "rle")) {
        # FIXME Need to rework this using `expect_type` instead.
        expect_is(
            object = counts(object, normalized = normalized),
            class = "matrix"
        )
        expect_null(assays(object)[[normalized]])
    }
})

test_that("Skipped DESeq2 transforms", {
    for (normalized in c("rlog", "vst")) {
        assays(object)[[normalized]] <- NULL
        expect_error(
            object = counts(object, normalized = normalized),
            regexp = normalized
        )
    }
})
