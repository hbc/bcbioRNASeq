context("plotPCACovariates")

skip_if_not_installed("DEGreport")

test_that("plotPCACovariates", {
    x <- plotPCACovariates(object)
    expect_is(x, "list")

    if (packageVersion("DEGreport") < "1.18") {
        expected <- c(
            "significantCovars",
            "plot",
            "corMatrix",
            "pcsMatrix",
            "scatterPlot",
            "effectsSignificantCovars"
        )
    } else {
        expected <- c(
            "plot",
            "corMatrix",
            "pcsMatrix",
            "scatterPlot",
            "significants"
        )
    }
    expect_identical(names(x), expected)
})
