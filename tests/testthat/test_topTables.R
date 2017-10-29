context("topTables")

test_that("topTables", {
    resTbl <- resultsTables(
        res,
        quiet = TRUE,
        summary = FALSE,
        write = FALSE)
    # Capture the knitr table output
    output <- capture_output(topTables(resTbl)) %>%
        strsplit("\\n") %>%
        .[[1]]
    # Expect 12 lines of knitr output
    expect_equal(length(output), 12)
    # Check for ensgene header
    expect_true(grepl("^\\|ensgene", output[[3]]))
})
