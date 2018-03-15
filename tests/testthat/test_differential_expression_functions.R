context("Differential Expression Functions")

# alphaSummary =================================================================
test_that("alphaSummary : bcbioRNASeq", {
    expect_warning(
        alphaSummary(bcb_small),
        "Empty design formula detected"
    )
    expect_is(suppressWarnings(alphaSummary(bcb_small)), "knitr_kable")
})

test_that("alphaSummary : DESeqDataSet", {
    x <- alphaSummary(dds_small)
    expect_is(x, "knitr_kable")
    expect_true(grepl("1e-06", x[[1L]]))
    # Upregulated genes
    expect_identical(
        x[[3L]],
        paste(
            "",
            "LFC > 0 (up)   ",
            "131, 26%          ",
            "112, 22%          ",
            "77, 15%           ",
            "40, 8%             ",
            "21, 4.2%          ",
            "",
            sep = "|"
        )
    )
    # Downregulated genes
    expect_identical(
        x[[4L]],
        paste(
            "",
            "LFC < 0 (down) ",
            "146, 29%          ",
            "128, 26%          ",
            "104, 21%          ",
            "67, 13%            ",
            "16, 3.2%          ",
            "",
            sep = "|"
        )
    )
})
