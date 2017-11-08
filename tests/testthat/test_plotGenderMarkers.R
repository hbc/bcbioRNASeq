context("plotGenderMarkers")

test_that("plotGenderMarkers", {
    # Use F1000 workflow example data (GSE65267) for unit tests
    # https://github.com/hbc/bcbioRNASeq/tree/f1000v1
    loadRemoteData(
        file.path(
            "https://github.com", "hbc", "bcbioRNASeq",
            "raw", "f1000v1", "data", "bcb.rda"),
        quiet = TRUE)
    p <- plotGenderMarkers(bcb)
    expect_is(p, "ggplot")
    rm(bcb)
})
