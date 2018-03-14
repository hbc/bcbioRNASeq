context("plotMeanSD")

test_that("plotGenderMarkers", {
    p <- plotMeanSD(bcb_small)
    expect_is(p, "ggplot")
})
