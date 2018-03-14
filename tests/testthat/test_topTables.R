context("topTables")

test_that("topTables", {
    resTbl <- resultsTables(
        res_small,
        summary = FALSE,
        write = FALSE
    )
    # Capture the knitr table output
    output <- capture_output(topTables(resTbl)) %>%
        strsplit("\\n") %>%
        .[[1L]]
    # Check for geneID column in header
    expect_true(grepl("^\\|geneID", output[[3L]]))
    # Check the output length
    expect_identical(length(output), 12L)
})
