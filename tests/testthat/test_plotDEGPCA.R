context("plotDEGPCA")

test_that("DESeqResults, DESeqTransform", {
    p <- plotDEGPCA(res_small, counts = rld_small)
    expect_is(p, "ggplot")
    expect_message(
        plotDEGPCA(res_small, counts = rld_small),
        "Plotting PCA using 4 genes"
    )
})
