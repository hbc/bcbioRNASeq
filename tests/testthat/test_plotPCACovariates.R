context("plotPCACovariates")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("Default", {
    p <- suppressWarnings(suppressMessages(
        plotPCACovariates(bcb)
    ))
    expect_is(p, "list")
    expect_identical(
        names(p),
        c("significantCovars",
          "plot",
          "corMatrix",
          "pcsMatrix",
          "scatterPlot",
          "effectsSignificantCovars")
    )
    # Check significant covariates
    expect_identical(
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
    expect_identical(
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

test_that("Defined metrics", {
    p <- plotPCACovariates(
        bcb,
        metrics = c("exonicRate", "intronicRate"))
    expect_identical(
        as.character(p[["significantCovars"]]),
        c("exonicRate", "intronicRate")
    )
    expect_identical(
        round(p[["effectsSignificantCovars"]], digits = 3L),
        c(exonicRate = 0.610,
          intronicRate = 0.610)
    )
})

test_that("Invalid parameters", {
    # Error on invalid column
    expect_error(
        plotPCACovariates(bcb, metrics = c("FOO", "BAR")),
        paste(
            "is_subset :",
            "The elements 'FOO', 'BAR' in col are not in",
            "colnames\\(metadata\\)."
        )
    )
    # More than 1 metric is required
    expect_error(
        plotPCACovariates(bcb, metrics = "exonicRate"),
        paste(
            "is_greater_than : length\\(col\\) are not all greater than 1L."
        )
    )
})

test_that("transformationLimit", {
    skip <- bcb
    assays(skip)[["rlog"]] <- NULL
    expect_warning(
        suppressMessages(
            plotPCACovariates(skip, normalized = "rlog")
        ),
        paste(
            "rlog counts not defined.",
            "Calculating and using log2 tmm counts on the fly instead."
        )
    )
    p <- suppressWarnings(suppressMessages(
        plotPCA(skip, normalized = "rlog")
    ))
    expect_is(p, "ggplot")
})
