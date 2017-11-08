context("txi")

test_that("txi", {
    txi <- txi(bcb)
    expect_equal(
        lapply(txi, class),
        list(abundance = "matrix",
             counts = "matrix",
             length = "matrix",
             countsFromAbundance = "character")
    )
})
