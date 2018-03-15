context("Gene Expression Functions")

genes <- head(rownames(bcb_small), 4L)
gene2symbol <- gene2symbol(bcb_small)



# plotGenderMarkers ============================================================
test_that("plotGenderMarkers : bcbioRNASeq", {
    p <- plotGenderMarkers(bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotGenderMarkers : DESeqDataSet", {
    p <- plotGenderMarkers(dds_small, interestingGroups = "treatment")
    expect_is(p, "ggplot")
})

test_that("plotGenderMarkers : DESeqTransform", {
    p <- plotGenderMarkers(rld_small, interestingGroups = "treatment")
    expect_is(p, "ggplot")
})



# plotGene =====================================================================
test_that("bcbioRNASeq", {
    p <- plotGene(bcb_small, genes = genes)
    expect_is(p, "ggplot")
})

test_that("List return", {
    list <- plotGene(bcb_small, genes = genes, return = "list")
    expect_is(list, "list")
    expect_true(
        lapply(list, function(x) is(x, "ggplot")) %>%
            unlist() %>%
            all()
    )
})

test_that("Markdown return", {
    gene <- gene2symbol[1L, "geneName", drop = TRUE]
    output <- capture.output(
        plotGene(bcb_small, genes = genes, return = "markdown")
    )
    expect_identical(
        output[[3L]],
        paste("##", gene)
    )
})

test_that("DESeqDataSet", {
    p <- plotGene(dds_small, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")
})

test_that("DESeqTransform", {
    p <- plotGene(rld_small, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")
})



# plotHeatmap ==================================================================
test_that("plotHeatmap : bcbioRNASeq", {
    genes <- head(rownames(bcb_small), n = 20L)
    p <- plotHeatmap(bcb_small, genes = genes)
    expect_is(p, "list")
    expect_identical(names(p), plotlist)
})

test_that("plotHeatmap : DESeqDataSet", {
    genes <- head(rownames(dds_small), n = 20L)
    p <- plotHeatmap(dds_small, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "list")
    expect_identical(names(p), plotlist)
})
