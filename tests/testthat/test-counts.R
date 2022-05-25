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
            expect_type(
                object = counts(object, normalized = normalized),
                type = "double"
            )
            expect_identical(
                object = counts(object, normalized = normalized),
                expected = assays(object)[[assay]]
            )
        }
    )
})

test_that("On the fly assay calculations", {
    for (normalized in c("tmm", "rle")) {
        expect_type(
            object = counts(object, normalized = normalized),
            type = "double"
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
