context("plotVolcano")

load(system.file(
    file.path("inst", "extdata", "res.rda"),
    package = "bcbioRNASeq"))

test_that("DESeqResults", {
    plot <- plotVolcano(res)
    expect_equal(
        class(plot),
        c("gg", "ggplot")
    )
})

test_that("data.frame", {
    df <- as.data.frame(res)
    plot <- plotVolcano(df)
    expect_equal(
        class(plot),
        c("gg", "ggplot")
    )
})
