context("plotGenderMarkers")

# Use F1000 workflow example data (GSE65267) for unit tests
# https://github.com/hbc/bcbioRNASeq/tree/f1000v1
loadRemoteData(
    file.path("http://bcbiornaseq.seq.cloud",
              "f1000v1",
              "data",
              "bcb.rda"),
    quiet = TRUE)

test_that("plotGenderMarkers", {
    p <- plotGenderMarkers(bcb)
    expect_is(p, "ggplot")
})
