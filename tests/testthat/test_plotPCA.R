context("plotPCA")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("plotPCA", {
    p <- plotPCA(bcb)
    expect_is(p, "ggplot")
})

test_that("transformationLimit", {
    skip <- bcb
    assays(skip)[["rlog"]] <- NULL
    expect_warning(
        plotPCA(skip, normalized = "rlog"),
        "rlog counts not defined. Using log2 tmm counts instead."
    )
    p <- suppressWarnings(
        plotPCA(skip, normalized = "rlog")
    )
    expect_is(p, "ggplot")
})
