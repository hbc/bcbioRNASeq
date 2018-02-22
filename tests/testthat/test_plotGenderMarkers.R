context("plotGenderMarkers")

url <- "http://bcbiornaseq.seq.cloud/f1000v1"
loadRemoteData(c(
    paste(url, "bcb.rda", sep = "/"),
    paste(url, "rld.rda", sep = "/"),
    path(url, "vst.rda", sep = "/")
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
