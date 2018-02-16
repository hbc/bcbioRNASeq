context("colData")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("bcbioRNASeq", {
    fix <- colData(bcb)
    fix[["batch"]] <- c(1L, 2L, 1L, 2L)
    colData(bcb) <- fix
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
