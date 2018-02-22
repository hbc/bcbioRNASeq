context("plotDEGHeatmap")

load(system.file("extdata/res.rda", package = "bcbioRNASeq"))
load(system.file("extdata/rld.rda", package = "bcbioRNASeq"))

test_that("plotDEGHeatmap", {
    p <- plotDEGHeatmap(object = res, counts = rld)
    expect_is(p, "list")
    expect_equal(
        names(p),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
})
