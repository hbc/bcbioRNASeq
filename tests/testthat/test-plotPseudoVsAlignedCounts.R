context("plotPseudoVsAlignedCounts")

test_that("bcbioRNASeq", {
    ## Correlation heatmap.
    p <- plotPseudoVsAlignedCounts(bcb)
    expect_s3_class(p, "pheatmap")
    ## Individual genes.
    ## Checking the most expressed aligned genes here.
    assay <- assays(bcb)[["aligned"]]
    genes <- names(tail(sort(rowSums(assay)), n = 2L))
    p <- plotPseudoVsAlignedCounts(bcb, genes = genes)
    expect_s3_class(p, "ggplot")
})
