context("plotCorrelationHeatmap")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("plotCorrelationHeatmap", {
    # Pearson (default)
    p <- plotCorrelationHeatmap(bcb)
    expect_equal(
        names(p),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
    # Spearman
    p <- plotCorrelationHeatmap(bcb, method = "spearman")
    expect_equal(
        names(p),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
    # Bad method
    expect_error(
        plotCorrelationHeatmap(bcb, method = "XXX"),
        "Supported methods: pearson, spearman"
    )
})

test_that("transformationLimit", {
    skip <- bcb
    assays(skip)[["rlog"]] <- NULL
    expect_warning(
        plotCorrelationHeatmap(skip, normalized = "rlog"),
        "rlog counts not defined. Using log2 tmm counts instead."
    )
    p <- suppressWarnings(
        plotCorrelationHeatmap(skip, normalized = "rlog")
    )
    expect_is(p, "list")
})
