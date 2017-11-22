context("plotGeneHeatmap")

bcb <- examples[["bcb"]]

test_that("plotGeneHeatmap", {
    genes <- counts(bcb)[1:20, ] %>% rownames()
    p <- plotGeneHeatmap(bcb, genes = genes)
    expect_equal(
        names(p),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
})
