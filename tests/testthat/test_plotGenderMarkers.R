context("plotGenderMarkers")

test_that("plotGenderMarkers", {
    p <- plotGenderMarkers(gse65267)
    expect_is(p, "ggplot")
})
