context("plotMA")

load(system.file("extdata/res.rda", package = "bcbioRNASeq"))
genes <- head(rownames(res))

test_that("DESeqResults", {
    p <- plotMA(res)
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
    p <- plotMA(res, genes = genes)
    expect_is(p, "ggplot")
})
