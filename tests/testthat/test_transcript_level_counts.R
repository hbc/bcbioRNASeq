context("Transcript-Level Counts")

bcb <- bcbioRNASeq(
    uploadDir = uploadDir,
    level = "transcripts",
    organism = "Mus musculus",
    ensemblRelease = 87L
)

test_that("counts", {
    # Valid: FALSE, tpm
    expect_is(counts(bcb, normalized = FALSE), "matrix")
    expect_is(counts(bcb, normalized = "tpm"), "matrix")
    # All other options are invalid.
    expect_error(counts(bcb, normalized = TRUE))
    expect_error(counts(bcb, normalized = "vst"))
})

test_that("Gene-level specific", {
    expect_error(
        object = plotDispEsts(bcb),
        regexp = "Gene-level counts are required"
    )
})

test_that("Show method", {
    object <- capture.output(show(bcb))
    expect_true(any(grepl("transcripts", object)))
})
