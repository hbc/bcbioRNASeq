context("plotGenderMarkers")

dir <- "http://bcbiornaseq.seq.cloud/f1000v1"
loadRemoteData(
    c(
        file.path(dir, "bcb.rda"),
        file.path(dir, "rld.rda"),
        file.path(dir, "vst.rda")
    ),
    quiet = TRUE)
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
