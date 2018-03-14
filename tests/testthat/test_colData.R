context("colData")

test_that("bcbioRNASeq", {
    data <- colData(bcb_small)
    data[["batch"]] <- seq_len(nrow(data))
    colData(bcb_small) <- data
    value <- colData(bcb_small)
    expect_true("batch" %in% colnames(value))
    expect_true(
        all(vapply(
            X = value,
            FUN = is.factor,
            FUN.VALUE = logical(1L)
        ))
    )
})
