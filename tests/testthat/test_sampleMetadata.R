context("sampleMetadata")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
load(system.file("extdata/dds.rda", package = "bcbioRNASeq"))
load(system.file("extdata/rld.rda", package = "bcbioRNASeq"))

test_that("bcbioRNASeq", {
    expect_true(all(lapply(sampleMetadata(bcb), is.factor)))
    expect_identical(
        sampleMetadata(bcb),
        as.data.frame(colData(bcb))
    )
})

test_that("DESeqDataSet", {
    expect_true(all(lapply(sampleMetadata(dds), is.factor)))
})

test_that("DESeqTransform", {
    expect_true(all(lapply(sampleMetadata(rld), is.factor)))
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
