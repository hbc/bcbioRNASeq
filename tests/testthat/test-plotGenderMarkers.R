context("plotGenderMarkers")

test_that("plotGenderMarkers", {
    expect_s3_class(
        object = plotGenderMarkers(object),
        class = "ggplot"
    )
})
