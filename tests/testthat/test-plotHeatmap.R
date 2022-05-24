test_that("plotHeatmap", {
    object <- nonzeroRowsAndCols(object)
    p <- plotHeatmap(object, scale = "row")
    expect_s3_class(p, "pheatmap")
})
