context("plotVolcano")

# TODO Define as global
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
