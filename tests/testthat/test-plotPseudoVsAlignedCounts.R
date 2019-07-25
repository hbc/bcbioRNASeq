context("plotPseudoVsAlignedCounts")

test_that("bcbioRNASeq", {
    ## Correlation heatmap.
    p <- plotPseudoVsAlignedCounts(bcb)
    expect_s3_class(p, "pheatmap")

    ## Individual genes.
    ## Checking the most expressed aligned genes here.
    genes <- assays(bcb) %>%
        .[["aligned"]] %>%
        rowSums() %>%
        sort() %>%
        tail(n = 2) %>%
        names()
    p <- plotPseudoVsAlignedCounts(bcb, genes = genes)
    expect_s3_class(p, "ggplot")
})
