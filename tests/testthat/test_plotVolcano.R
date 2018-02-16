context("plotVolcano")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))
load(system.file(
    file.path("extdata", "res.rda"),
    package = "bcbioRNASeq"))

gene2symbol <- gene2symbol(bcb)

test_that("DESeqResults", {
    p <- plotVolcano(res, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")

    # Label the top genes
    p <- plotVolcano(res, ntop = 5L, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")

    # Label specific genes
    genes <- rownames(res) %>% head()
    p <- plotVolcano(res, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")
})

test_that("Data frame", {
    df <- as.data.frame(res)
    p <- plotVolcano(df)
    expect_is(p, "ggplot")
})
