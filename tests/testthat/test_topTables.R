context("topTables")

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
    expect_identical(length(output), 12L)
})
