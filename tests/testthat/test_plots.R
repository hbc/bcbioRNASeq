test_that("Plots: MA plot", {
    p <- plot_ma(res)

    # Check geom classes
    geom_type <- sapply(p[["layers"]], function(x) class(x[["geom"]])[[1L]])
    expect_identical(geom_type, c("GeomPoint", "GeomLogticks"))

    # Check plot labels
    expect_identical(p[["labels"]][["y"]],
                     "log2 fold change")
    expect_identical(p[["labels"]][["x"]],
                     "mean expression across all samples")
})
