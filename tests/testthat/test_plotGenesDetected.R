context("plotGenesDetected")

test_that("plotGenesDetected", {
    p <- plotGenesDetected(bcb_small)
    expect_is(p, "ggplot")
})
