context("plotMeanSD")

test_that("plotGenderMarkers", {
    p <- plotMeanSD(bcb)
    expect_true(is(p, "ggplot"))
})
