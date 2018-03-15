context("Deprecated Functions")



# sampleMetadata ===============================================================
test_that("bcbioRNASeq", {
    x <- suppressWarnings(sampleMetadata(bcb_small))
    expect_true(all(vapply(
        X = x,
        FUN = is.factor,
        FUN.VALUE = logical(1L)
    )))
    expect_identical(x, colData(bcb_small))
})

test_that("DESeqDataSet", {
    x <- suppressWarnings(sampleMetadata(dds_small))
    expect_true(all(vapply(
        X = x,
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
    bcb <- bcb_small
    x <- suppressWarnings(sampleMetadata(bcb))
    x[["batch"]] <- seq_len(nrow(x))
    suppressWarnings(sampleMetadata(bcb) <- x)
    expect_s4_class(bcb, "bcbioRNASeq")
    value <- suppressWarnings(sampleMetadata(bcb))
    expect_true("batch" %in% colnames(value))
    expect_true(
        all(vapply(
            X = value,
            FUN = is.factor,
            FUN.VALUE = logical(1L)
        ))
    )
})
