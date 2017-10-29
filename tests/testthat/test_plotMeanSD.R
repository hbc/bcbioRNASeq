context("plotMeanSD")

test_that("plotGenderMarkers", {
    p <- plotMeanSD(bcb)
    expect_is(p, "ggplot")
})
