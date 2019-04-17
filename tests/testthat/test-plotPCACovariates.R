context("plotPCACovariates")

test_that("plotPCACovariates", {
    x <- plotPCACovariates(object)
    expect_is(x, "list")
    expect_identical(
        names(x),
        c(
            "plot",
            "corMatrix",
            "pcsMatrix",
            "scatterPlot",
            "significants"
        )
    )
})
