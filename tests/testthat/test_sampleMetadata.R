context("sampleMetadata")

test_that("bcbioRNASeq", {
    expect_true(all(vapply(
        X = sampleMetadata(bcb_small),
        FUN = is.factor,
        FUN.VALUE = logical(1L)
    )))
    expect_identical(
        sampleMetadata(bcb_small),
        as.data.frame(colData(bcb_small))
    )
})

test_that("DESeqDataSet", {
    expect_true(all(vapply(
        X = sampleMetadata(dds_small),
        FUN = is.factor,
        FUN.VALUE = logical(1L)
    )))
})

test_that("DESeqTransform", {
    expect_true(all(vapply(
        X = sampleMetadata(rld_small),
        FUN = is.factor,
        FUN.VALUE = logical(1L)
    )))
})

test_that("Assignment method", {
    data <- sampleMetadata(bcb_small)
    data[["batch"]] <- c(1L, 2L, 1L, 2L)
    sampleMetadata(bcb_small) <- data
    value <- sampleMetadata(bcb_small)
    expect_identical(dim(value), c(4L, 5L))
    expect_true(
        all(vapply(
            X = value,
            FUN = is.factor,
            FUN.VALUE = logical(1L)
        ))
    )
})
