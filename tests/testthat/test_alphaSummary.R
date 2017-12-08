context("alphaSummary")

load(system.file(
    file.path("inst", "extdata", "bcb.rda"),
    package = "bcbioRNASeq"))
load(system.file(
    file.path("inst", "extdata", "dds.rda"),
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
               "170, 56%         |",
               "166, 55%         |",
               "157, 52%         |",
               "152, 50%         |",
               "136, 45%         |")
    )
    # Downregulated genes
    expect_identical(
        summary[[4L]],
        paste0("|LFC < 0 (down) |",
               "0, 0%            |",
               "0, 0%            |",
               "0, 0%            |",
               "0, 0%            |",
               "0, 0%            |")
    )
})
