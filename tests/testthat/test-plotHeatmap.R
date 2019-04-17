context("plotHeatmap")

test_that("plotHeatmap", {
    genes <- head(rownames(object), n = 100L)
    p <- plotHeatmap(object[genes, ])
    expect_s3_class(p, "pheatmap")
})
