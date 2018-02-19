context("plotRRNAMappingRate")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("bcbioRNASeq", {
    p <- plotRRNAMappingRate(bcb)
    expect_is(p, "ggplot")
})

test_that("Data frame", {
    data <- metrics(bcb)
    p <- plotRRNAMappingRate(bcb)
    expect_is(p, "ggplot")
})

test_that("Legacy rRnaRate column", {
    data <- metrics(bcb)
    data[["rRnaRate"]] <- data[["rrnaRate"]]
    data[["rrnaRate"]] <- NULL
    expect_warning(
        plotRRNAMappingRate(data),
        "`rrnaRate` is missing from `metrics\\(\\)`"
    )
})
