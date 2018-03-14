context("selectSamples")

# TODO Define this as a global
dim <- c(505L, 2L)

test_that("bcbioRNASeq", {
    x <- selectSamples(bcb_small, group = "ko", transform = TRUE)
    expect_identical(dim(x), dim)
    x <- selectSamples(bcb_small, group = "ko", transform = FALSE)
    expect_identical(dim(x), dim)
})

test_that("DESeqDataSet", {
    x <- selectSamples(dds_small, group = "ko")
    expect_identical(dim(x), dim)
})
