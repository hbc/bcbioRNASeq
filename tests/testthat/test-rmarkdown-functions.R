context("R Markdown Functions")



# topTables ====================================================================
test_that("topTables : resultsTables list", {
    x <- resultsTables(
        results = res_small,
        counts = dds_small,
        summary = FALSE,
        write = FALSE
    )
    # Capture the knitr table output
    output <- capture_output(topTables(x)) %>%
        strsplit("\\n") %>%
        .[[1L]]
    # Check for geneID column in header
    expect_true(grepl("geneID", output[[3L]]))

    # Coding mode
    output <- capture_output(topTables(x, coding = TRUE))
    expect_true(grepl("protein_coding", output))
})

test_that("topTables : DESeqResults", {
    output <- capture.output(topTables(res_small))
    expect_true(grepl("^\\|geneID", output[[3L]]))
})
