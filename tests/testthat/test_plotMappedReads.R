context("plotMappedReads")

test_that("plotMappedReads", {
    p <- plotMappedReads(bcb)
    expect_true(is(p, "ggplot"))
})
