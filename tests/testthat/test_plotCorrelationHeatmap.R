context("plotCorrelationHeatmap")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("plotCorrelationHeatmap", {
    # Pearson (default)
    p <- plotCorrelationHeatmap(bcb)
    expect_identical(
        names(p),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
    # Spearman
    p <- plotCorrelationHeatmap(bcb, method = "spearman")
    expect_identical(
        names(p),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
    # Bad method
    expect_error(
        plotCorrelationHeatmap(bcb, method = "XXX"),
        paste(
            "is_subset :",
            "The element 'XXX' in method is not in",
            "c\\(\"pearson\", \"spearman\"\\)."
        )
    )
})

test_that("transformationLimit", {
    skip <- bcb
    assays(skip)[["rlog"]] <- NULL
    expect_warning(
        plotCorrelationHeatmap(skip, normalized = "rlog"),
        paste(
            "rlog counts not defined.",
            "Calculating and using log2 tmm counts on the fly instead."
        )
    )
    p <- suppressWarnings(plotCorrelationHeatmap(skip, normalized = "rlog"))
    expect_is(p, "list")
})
