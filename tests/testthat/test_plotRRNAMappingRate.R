context("plotRRNAMappingRate")

test_that("bcbioRNASeq", {
    p <- plotRRNAMappingRate(bcb_small)
    expect_is(p, "ggplot")
})

test_that("Data frame", {
    data <- metrics(bcb_small)
    p <- plotRRNAMappingRate(bcb_small)
    expect_is(p, "ggplot")
})

test_that("Legacy rRnaRate column", {
    data <- metadata(bcb_small)[["metrics"]]
    data[["rRnaRate"]] <- data[["rrnaRate"]]
    data[["rrnaRate"]] <- NULL
    metadata(bcb_small)[["metrics"]] <- data
    expect_warning(
        plotRRNAMappingRate(bcb_small),
        "`rrnaRate` is missing from `metrics\\(\\)`"
    )
})
