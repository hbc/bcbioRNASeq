context("plotPCACovariates")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("default", {
    p <- suppressWarnings(suppressMessages(
        plotPCACovariates(bcb)
    ))
    expect_is(p, "list")
    expect_equal(
        names(p),
        c("significantCovars",
          "plot",
          "corMatrix",
          "pcsMatrix",
          "scatterPlot",
          "effectsSignificantCovars")
    )
    # Check significant covariates
    expect_equal(
        p[["significantCovars"]] %>%
            na.omit() %>%
            as.character() %>%
            sort(),
        c("duplicationRateOfMapped",
          "exonicRate",
          "intergenicRate",
          "intronicRate",
          "rrna")
    )
    expect_equal(
        p[["effectsSignificantCovars"]] %>%
            na.omit() %>%
            .[. > 0L] %>%
            .[(sort(names(.)))] %>%
            round(digits = 3L),
        c(duplicationRateOfMapped = 0.610,
          exonicRate = 0.610,
          intergenicRate = 0.610,
          intronicRate = 0.610,
          rrna = 0.279)
    )
})

test_that("defined metrics", {
    p <- plotPCACovariates(
        bcb,
        metrics = c("exonicRate", "intronicRate"))
    expect_equal(
        as.character(p[["significantCovars"]]),
        c("exonicRate", "intronicRate")
    )
    expect_equal(
        round(p[["effectsSignificantCovars"]], digits = 3L),
        c(exonicRate = 0.610,
          intronicRate = 0.610)
    )
})

test_that("invalid parameters", {
    # Error on invalid column
    expect_error(
        plotPCACovariates(bcb, metrics = c("FOO", "BAR")),
        "Failed to select valid `metrics` columns to plot"
    )
    # More than 1 metric is required
    expect_error(
        plotPCACovariates(bcb, metrics = "exonicRate"),
        "`degCovariates\\(\\)` requires at least 2 metadata columns"
    )
})

test_that("transformationLimit", {
    skip <- bcb
    assays(skip)[["rlog"]] <- NULL
    expect_warning(
        suppressMessages(
            plotPCACovariates(skip, normalized = "rlog")
        ),
        "rlog counts not defined. Using log2 tmm counts instead."
    )
    p <- suppressWarnings(suppressMessages(
        plotPCA(skip, normalized = "rlog")
    ))
    expect_is(p, "ggplot")
})
