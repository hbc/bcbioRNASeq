context("alphaSummary")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))
load(system.file(
    file.path("extdata", "dds.rda"),
    package = "bcbioRNASeq"))

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
    expect_true(grepl("1e-06", summary[[1L]]))
    # Upregulated genes
    expect_identical(
        summary[[3L]],
        paste0("|LFC > 0 (up)   |",
               "2, 0.66%         |",
               "0, 0%            |",
               "0, 0%            |",
               "0, 0%            |",
               "0, 0%            |")
    )
    # Downregulated genes
    expect_identical(
        summary[[4L]],
        paste0("|LFC < 0 (down) |",
               "2, 0.66%         |",
               "0, 0%            |",
               "0, 0%            |",
               "0, 0%            |",
               "0, 0%            |")
    )
})
