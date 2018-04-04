context("R Markdown Functions")



# prepareRNASeqTemplate ========================================================
test_that("prepareRNASeqTemplate", {
    files <- c(
        "_footer.Rmd",
        "_header.Rmd",
        "_output.yaml",
        "_setup.R",
        "bibliography.bib"
    )
    expect_silent(prepareRNASeqTemplate())
    expect_true(all(file.exists(files)))
    unlink(files)
})



# topTables ====================================================================
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
    expect_identical(length(output), 108L)
})
