context("alphaSummary")

test_that("bcbioRNASeq", {
    expect_warning(
        alphaSummary(bcb),
        "Empty DESeqDataSet design formula detected"
    )
    expect_is(suppressWarnings(alphaSummary(bcb)), "knitr_kable")
})

test_that("DESeqData", {
    summary <- alphaSummary(dds)
    expect_is(summary, "knitr_kable")
    expect_true(grepl("1e-06", summary[[1]]))
    expect_true(grepl("\\|2, 0.66%", summary[[3]]))
})
