context("deprecated")

bcb <- examples[["bcb"]]

test_that("download", {
    expect_warning(
        download(),
        "Use 'prepareRNASeqTemplate' instead."
    )
})

test_that("loadRNASeqRun", {
    uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
    expect_warning(
        loadRNASeqRun(uploadDir),
        "Use 'loadRNASeq' instead."
    )
})

test_that("plotGeneDetectionSaturation", {
    expect_warning(
        plotGeneDetectionSaturation(bcb),
        "Use 'plotGeneSaturation' instead."
    )
})

test_that("plotDispersion", {
    expect_warning(
        plotDispersion(bcb),
        "Use 'plotDispEsts' instead."
    )
})
