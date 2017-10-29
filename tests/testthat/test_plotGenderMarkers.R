context("plotGenderMarkers")

test_that("plotGenderMarkers", {
    file.path(
        "https://github.com",
        "hbc",
        "bcbioRNASeq",
        "raw",
        "f1000v1",
        "data",
        "bcb.rda") %>%
        loadRemoteData(quiet = TRUE)
    p <- plotGenderMarkers(bcb)
    expect_true(is(p, "ggplot"))
})
