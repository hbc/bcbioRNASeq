context("plotGenesDetected")

test_that("plotGenesDetected", {
    p <- plotGenesDetected(bcb)
    expect_is(p, "ggplot")
})
