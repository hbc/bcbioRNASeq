context("plotGenesDetected")

test_that("plotGenesDetected", {
    p <- plotGenesDetected(bcb)
    expect_true(is(p, "ggplot"))
})
