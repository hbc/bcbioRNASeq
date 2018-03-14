context("plotGene")

# TODO Define these as globals
genes <- rownames(bcb_small)[1L:4L]
gene2symbol <- gene2symbol(bcb_small)

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
    output <- capture.output(
        plotGene(bcb_small, genes = genes, return = "markdown")
    )
    expect_identical(output[[3L]], "## Rgs20")
})

test_that("DESeqDataSet", {
    p <- plotGene(dds_small, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")
})

test_that("DESeqTransform", {
    p <- plotGene(rld_small, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")
})
