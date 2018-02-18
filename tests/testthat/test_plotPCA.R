context("plotPCA")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("plotPCA", {
    p <- plotPCA(bcb)
    expect_is(p, "ggplot")
})

test_that("transformationLimit", {
    skip <- bcb
    assays(skip)[["rlog"]] <- NULL
    expect_warning(
        plotPCA(skip, normalized = "rlog"),
        paste(
            "rlog counts not defined.",
            "Calculating and using log2 tmm counts on the fly instead."
        )
    )
    p <- suppressWarnings(
        plotPCA(skip, normalized = "rlog")
    )
    expect_is(p, "ggplot")
})
