context("plotGenderMarkers")

# TODO Migrate to using this example dataset as the main unit test set
url <- paste(cacheURL, "f1000v1", sep = "/")
loadRemoteData(c(
    paste(url, "bcb.rda", sep = "/"),
    paste(url, "rld.rda", sep = "/"),
    paste(url, "vst.rda", sep = "/")
))
dds <- bcbio(bcb, "DESeqDataSet")

test_that("bcbioRNASeq", {
    p <- plotGenderMarkers(bcb)
    expect_is(p, "ggplot")
})

test_that("DESeqDataSet", {
    p <- plotGenderMarkers(dds, interestingGroups = "group")
    expect_is(p, "ggplot")
})

test_that("DESeqTransform", {
    p <- plotGenderMarkers(rld, interestingGroups = "group")
    expect_is(p, "ggplot")
    p <- plotGenderMarkers(vst, interestingGroups = "group")
    expect_is(p, "ggplot")
})
