## Can hit this cryptic ggplot2 error on macOS unless we call at the top.
## > no applicable method for 'depth' applied to an object of class "NULL"
## https://github.com/tidyverse/ggplot2/issues/2514
grid::current.viewport()

object <- bcb

test_that("method", {
    ## Pearson.
    expect_s3_class(
        object = plotCorrelationHeatmap(object, method = "pearson"),
        class = "pheatmap"
    )
    ## Spearman.
    expect_s3_class(
        object = plotCorrelationHeatmap(object, method = "spearman"),
        class = "pheatmap"
    )
    ## Bad method.
    expect_error(
        object = plotCorrelationHeatmap(object, method = "XXX"),
        regexp = "'arg' should be one of"
    )
})

test_that("Fast mode", {
    expect_error(
        object = plotCorrelationHeatmap(bcb_fast),
        regexp = "fast mode"
    )
})
