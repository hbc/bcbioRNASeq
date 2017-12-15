context("plotGenderMarkers")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("plotGenderMarkers", {
    # Current working example doesn't contain the dimorphic genes
    expect_warning(
        plotGenderMarkers(bcb),
        "Missing gender markers in count matrix"
    )
    p <- suppressWarnings(plotGenderMarkers(bcb))
    expect_is(p, "NULL")
})
