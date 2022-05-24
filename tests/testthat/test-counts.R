test_that("Slotted assays", {
    mapply(
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
        FUN = function(normalized, assay) {
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
        },
        SIMPLIFY = FALSE
    )
})

test_that("On the fly assay calculations", {
    for (normalized in c("tmm", "rle")) {
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
