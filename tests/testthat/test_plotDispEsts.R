context("plotDispEsts")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("plotDispEsts", {
    p <- plotDispEsts(bcb)
    expect_is(p, "list")
    expect_identical(
        names(p),
        c("rect", "text")
    )
})
