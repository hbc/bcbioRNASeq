context("plotDEGenePCA")

load(system.file(
    file.path("extdata", "res.rda"),
    package = "bcbioRNASeq"))
load(system.file(
    file.path("extdata", "rld.rda"),
    package = "bcbioRNASeq"))

test_that("DESeqResults, DESeqTransform", {
    p <- suppressMessages(
        plotDEGenePCA(res, counts = rld)
    )
    expect_is(p, "ggplot")
    expect_message(
        plotDEGenePCA(res, counts = rld),
        "Plotting PCA using 4 genes"
    )
})
