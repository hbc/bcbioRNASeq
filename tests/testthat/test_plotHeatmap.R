context("plotHeatmap")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("plotHeatmap", {
    genes <- counts(bcb)[1L:20L, ] %>%
        rownames()
    p <- plotHeatmap(bcb, genes = genes)
    expect_equal(
        names(p),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
})
