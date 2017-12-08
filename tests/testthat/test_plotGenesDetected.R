context("plotGenesDetected")

load(system.file(
    file.path("inst", "extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("plotGenesDetected", {
    p <- plotGenesDetected(bcb)
    expect_is(p, "ggplot")
})
