context("plotGenderMarkers")

test_that("bcbioRNASeq", {
    rownames(object)[seq_len(2L)] <-
        c("ENSMUSG00000086503", "ENSMUSG00000069045")
    p <- plotGenderMarkers(object)
    expect_s3_class(p, "ggplot")
})
