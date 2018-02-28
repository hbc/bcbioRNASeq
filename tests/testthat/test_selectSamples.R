context("selectSamples")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
load(system.file("extdata/dds.rda", package = "bcbioRNASeq"))

dim <- c(505L, 2L)

test_that("bcbioRNASeq", {
    x <- selectSamples(bcb, group = "ko", transform = TRUE)
    expect_identical(dim(x), dim)
    expect_identical(
        names(assays(x)),
        c("raw", "normalized", "tpm", "tmm", "rlog", "vst")
    )
    expect_identical(
        rownames(metrics(x)),
        c("group2_1", "group2_2")
    )

    x <- selectSamples(bcb, group = "ko", transform = FALSE)
    expect_identical(dim(x), dim)
    expect_identical(
        names(assays(x)),
        c("raw", "normalized", "tpm", "tmm")
    )
    expect_identical(
        rownames(metrics(x)),
        c("group2_1", "group2_2")
    )
})

test_that("DESeqDataSet", {
    x <- selectSamples(dds, group = "ko")
    expect_identical(dim(x), dim)
})
