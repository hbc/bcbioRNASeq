context("plotGenderMarkers")

test_that("bcbioRNASeq", {
    p <- plotGenderMarkers(bcb_small)
    expect_is(p, "ggplot")
})

test_that("DESeqDataSet", {
    p <- plotGenderMarkers(dds_small, interestingGroups = "group")
    expect_is(p, "ggplot")
})

test_that("DESeqTransform", {
    p <- plotGenderMarkers(rld_small, interestingGroups = "group")
    expect_is(p, "ggplot")
    p <- plotGenderMarkers(vst, interestingGroups = "group")
    expect_is(p, "ggplot")
})
