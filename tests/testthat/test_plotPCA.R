context("plotPCA")

load(system.file(
    file.path("inst", "extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("plotPCA", {
    p <- plotPCA(bcb)
    expect_is(p, "ggplot")
})
