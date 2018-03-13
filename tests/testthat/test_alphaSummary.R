context("alphaSummary")

test_that("bcbioRNASeq", {
    expect_warning(
        alphaSummary(bcb_small),
        "Internal DESeqDataSet has an empty design formula"
    )
    expect_is(suppressWarnings(alphaSummary(bcb)), "knitr_kable")
})

test_that("DESeqData", {
    x <- alphaSummary(dds)
    expect_is(x, "knitr_kable")
    expect_true(grepl("1e-06", x[[1L]]))
    # Upregulated genes
    expect_identical(
        x[[3L]],
        paste0(
            "|LFC > 0 (up)   |",
            "2, 0.66%         |",
            "0, 0%            |",
            "0, 0%            |",
            "0, 0%            |",
            "0, 0%            |")
    )
    # Downregulated genes
    expect_identical(
        x[[4L]],
        paste0(
            "|LFC < 0 (down) |",
            "2, 0.66%         |",
            "0, 0%            |",
            "0, 0%            |",
            "0, 0%            |",
            "0, 0%            |")
    )
})
