context("plotDEGHeatmap")

test_that("plotDEGHeatmap", {
    p <- plotDEGHeatmap(res_small, counts = rld_small)
    expect_is(p, "list")
    expect_identical(
        names(p),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
})
