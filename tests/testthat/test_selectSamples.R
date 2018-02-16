context("selectSamples")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))
load(system.file(
    file.path("extdata", "dds.rda"),
    package = "bcbioRNASeq"))

dim <- c(505L, 2L)

test_that("bcbioRNASeq", {
    x <- suppressMessages(
        selectSamples(bcb, group = "ko", transform = TRUE)
    )
    expect_identical(dim(x), dim)
    x <- suppressMessages(
        selectSamples(bcb, group = "ko", transform = FALSE)
    )
    expect_identical(dim(x), dim)
})

test_that("DESeqDataSet", {
    x <- selectSamples(dds, group = "ko")
    expect_identical(dim(x), dim)
})
