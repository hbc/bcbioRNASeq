context("Markdown")



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
