context("plotHeatmap")

bcb <- examples[["bcb"]]

test_that("plotHeatmap", {
    genes <- counts(bcb)[1:20, ] %>%
        rownames()
    p <- plotHeatmap(bcb, genes = genes)
    expect_equal(
        names(p),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
})
