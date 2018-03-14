context("plotHeatmap")

# TODO Define as globals
gene2symbol <- gene2symbol(bcb_small)

test_that("bcbioRNASeq", {
    genes <- head(rownames(bcb_small), n = 20L)
    p <- plotHeatmap(bcb_small, genes = genes)
    expect_is(p, "list")
    expect_identical(names(p), plotlist)
})

test_that("DESeqDataSet", {
    genes <- head(rownames(dds_small), n = 20L)
    p <- plotHeatmap(dds_small, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "list")
    expect_identical(names(p), plotlist)
})
