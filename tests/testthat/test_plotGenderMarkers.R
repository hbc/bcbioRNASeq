context("plotGenderMarkers")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

# TODO Improve working example to include dimorphic genes

test_that("Missing markers in minimal example", {
    # Current working example doesn't contain the dimorphic genes
    expect_error(
        plotGenderMarkers(bcb),
        paste(
            "is_subset :",
            "The elements 'ENSMUSG00000056673'"
        )
    )
})
