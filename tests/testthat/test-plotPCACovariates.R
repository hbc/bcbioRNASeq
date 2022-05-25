test_that("plotPCACovariates", {
    x <- plotPCACovariates(object)
    expect_type(x, "list")
    expected <- c(
        "plot",
        "corMatrix",
        "pcsMatrix",
        "scatterPlot",
        "significants"
    )
    expect_named(x, expected)
})
