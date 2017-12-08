context("plot53Bias")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("plot53Bias", {
    p <- plot53Bias(bcb)
    expect_is(p, "ggplot")
})
