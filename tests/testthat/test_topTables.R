context("topTables")

load(system.file("extdata/res.rda", package = "bcbioRNASeq"))

test_that("topTables", {
    resTbl <- resultsTables(
        res,
        summary = FALSE,
        write = FALSE)
    # Capture the knitr table output
    output <- capture_output(topTables(resTbl)) %>%
        strsplit("\\n") %>%
        .[[1L]]
    # Check for ensgene column in header
    expect_true(grepl("^\\|ensgene", output[[3L]]))
    # Check the output length
    expect_equal(length(output), 12L)
})
