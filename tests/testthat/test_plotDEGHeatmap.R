context("plotDEGHeatmap")

test_that("plotDEGHeatmap", {
    p <- plotDEGHeatmap(object = res, counts = rld, quiet = TRUE)
    expect_is(p, "list")
    expect_equal(
        names(p),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
})
