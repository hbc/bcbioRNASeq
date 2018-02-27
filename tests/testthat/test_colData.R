context("colData")

test_that("bcbioRNASeq", {
    data <- colData(bcb)
    data[["batch"]] <- c(1L, 2L, 1L, 2L)
    colData(bcb) <- data
    value <- colData(bcb)
    expect_identical(dim(value), c(4L, 5L))
    expect_true(
        all(vapply(
            X = value,
            FUN = is.factor,
            FUN.VALUE = logical(1L)
        ))
    )
})
