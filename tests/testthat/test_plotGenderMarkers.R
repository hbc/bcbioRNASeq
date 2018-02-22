context("plotGenderMarkers")

dir <- "http://bcbiornaseq.seq.cloud/f1000v1"
loadRemoteData(
    c(path(dir, "bcb.rda"), path(dir, "rld.rda"), path(dir, "vst.rda"))
)
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
