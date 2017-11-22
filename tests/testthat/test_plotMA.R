context("plotMA")

res <- examples[["res"]]

test_that("plotMA", {
    p <- plotMA(res)

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
