context("plotVolcano")

# TODO Define as global
gene2symbol <- gene2symbol(bcb_small)

test_that("DESeqResults", {
    p <- plotVolcano(res_small, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")

    # Label the top genes
    p <- plotVolcano(res_small, ntop = 5L, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")

    # Label specific genes
    genes <- rownames(res_small) %>% head()
    p <- plotVolcano(res_small, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")
})
