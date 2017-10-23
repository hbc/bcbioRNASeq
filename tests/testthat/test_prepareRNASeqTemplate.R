context("prepareRNASeqTemplate")

test_that("prepareRNASeqTemplate", {
    files <- c(
        "_footer.Rmd",
        "_header.Rmd",
        "_output.yaml",
        "bibliography.bib",
        "setup.R")
    expect_silent(prepareRNASeqTemplate())
    expect_true(all(file.exists(files)))
    unlink(files)
})
