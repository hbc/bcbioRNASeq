context("plotHeatmap")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
load(system.file("extdata/dds.rda", package = "bcbioRNASeq"))
load(system.file("extdata/rld.rda", package = "bcbioRNASeq"))

gene2symbol <- gene2symbol(bcb)
plotlist <- c("tree_row", "tree_col", "kmeans", "gtable")

test_that("bcbioRNASeq", {
    genes <- head(rownames(bcb), n = 20L)
    p <- plotHeatmap(bcb, genes = genes)
    expect_is(p, "list")
    expect_identical(names(p), plotlist)
})

test_that("DESeqDataSet", {
    genes <- head(rownames(dds), n = 20L)
    p <- plotHeatmap(dds, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "list")
    expect_identical(names(p), plotlist)
})
