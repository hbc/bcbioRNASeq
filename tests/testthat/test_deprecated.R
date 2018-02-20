context("deprecated")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("download", {
    expect_warning(
        download(),
        "Use 'prepareRNASeqTemplate' instead."
    )
})

test_that("loadRNASeqRun", {
    expect_warning(
        loadRNASeqRun(),
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

test_that("plotGeneHeatmap", {
    expect_warning(
        plotGeneHeatmap(bcb),
        "Use 'plotHeatmap' instead."
    )
})

test_that("txi", {
    expect_warning(
        txi(bcb),
        "Use 'bcbio\\(object, \"tximport\"\\)' instead."
    )
})
