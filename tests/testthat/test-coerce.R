context("coerce")

test_that("bcbioRNASeq to DESeqDataSet", {
    x <- as(object, "DESeqDataSet")
    expect_s4_class(x, "DESeqDataSet")
})

test_that("bcbioRNASeq to DESeqTransform", {
    x <- as(object, "DESeqTransform")
    expect_s4_class(x, "DESeqTransform")
})
