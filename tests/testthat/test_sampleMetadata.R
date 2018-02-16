context("sampleMetadata")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))
load(system.file(
    file.path("extdata", "dds.rda"),
    package = "bcbioRNASeq"))
load(system.file(
    file.path("extdata", "rld.rda"),
    package = "bcbioRNASeq"))

test_that("bcbioRNASeq", {
    expect_identical(
        as.data.frame(colData(bcb)),
        sampleMetadata(bcb)
    )
})

test_that("DESeqDataSet", {
    expect_identical(
        as.data.frame(colData(dds)),
        sampleMetadata(dds)
    )
})

test_that("DESeqTransform", {
    expect_identical(
        as.data.frame(colData(rld)),
        sampleMetadata(rld)
    )
})

test_that("Assignment method", {
    data <- sampleMetadata(bcb)
    data[["batch"]] <- c(1L, 2L, 1L, 2L)
    sampleMetadata(bcb) <- data
    value <- sampleMetadata(bcb)
    expect_identical(dim(value), c(4L, 5L))
    expect_true(
        all(vapply(
            X = value,
            FUN = is.factor,
            FUN.VALUE = logical(1L)
        ))
    )
})
