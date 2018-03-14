context("plotMA")

# TODO Define as a global
genes <- head(rownames(res_small))

test_that("DESeqResults", {
    p <- plotMA(res_small)
    expect_is(p, "ggplot")

    # Check geom classes
    geomtype <- sapply(p[["layers"]], function(x) class(x[["geom"]])[[1L]])
    expect_identical(
        geomtype, c("GeomPoint", "GeomLogticks"))

    # Check plot labels
    expect_identical(
        p[["labels"]][["y"]],
        "log2 fold change")
    expect_identical(
        p[["labels"]][["x"]],
        "mean expression across all samples")
})

test_that("Gene labels", {
    p <- plotMA(res_small, genes = genes)
    expect_is(p, "ggplot")
})
