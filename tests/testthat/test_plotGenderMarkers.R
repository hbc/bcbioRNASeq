context("plotGenderMarkers")

test_that("plotGenderMarkers", {
    p <- plotGenderMarkers(bcb)
    expect_true(is(p, "ggplot"))
})
