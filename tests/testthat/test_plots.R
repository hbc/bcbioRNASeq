test_that("Plots: MA plot", {
    p <- plot_ma(res)

    # Check geom classes
    geom_type <- sapply(p[["layers"]], function(x) class(x[["geom"]])[[1L]])
    expect_identical(geom_type, c("GeomPoint", "GeomLogticks"))

    # Check plot labels
    expect_identical(p[["labels"]][["y"]],
                     expression(log[2] * " fold change"))  # nolint
    expect_identical(p[["labels"]][["x"]],
                     "mean expression across all samples")
})
