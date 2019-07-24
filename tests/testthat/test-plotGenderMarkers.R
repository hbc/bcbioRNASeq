context("plotGenderMarkers")

## Current minimal example doesn't contain dimorphic genes.
test_that("bcbioRNASeq", {
    expect_message(
        object = plotGenderMarkers(object),
        regexp = "ENSMUSG"
    )
})
