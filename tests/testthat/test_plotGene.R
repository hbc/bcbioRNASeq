context("plotGene")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
load(system.file("extdata/dds.rda", package = "bcbioRNASeq"))
load(system.file("extdata/rld.rda", package = "bcbioRNASeq"))

genes <- rownames(bcb)[1L:4L]
gene2symbol <- gene2symbol(bcb)

test_that("bcbioRNASeq", {
    p <- plotGene(bcb, genes = genes)
    expect_is(p, "ggplot")
})

test_that("List return", {
    list <- plotGene(bcb, genes = genes, return = "list")
    expect_is(list, "list")
    expect_true(
        lapply(list, function(x) is(x, "ggplot")) %>%
            unlist() %>%
            all()
    )
})

test_that("Markdown return", {
    output <- capture.output(
        plotGene(bcb, genes = genes, return = "markdown")
    )
    expect_identical(output[[3L]], "## Rgs20")
})

test_that("DESeqDataSet", {
    p <- plotGene(dds, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")
})

test_that("DESeqTransform", {
    p <- plotGene(rld, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")
})
