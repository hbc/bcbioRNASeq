context("plotDEGPCA")

load(system.file("extdata/res.rda", package = "bcbioRNASeq"))
load(system.file("extdata/rld.rda", package = "bcbioRNASeq"))

test_that("DESeqResults, DESeqTransform", {
    p <- plotDEGPCA(res, counts = rld)
    expect_is(p, "ggplot")
    expect_message(
        plotDEGPCA(res, counts = rld),
        "Plotting PCA using 4 genes"
    )
})
