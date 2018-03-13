context("design")

test_that("bcbioRNASeq", {
    # Not identical because of environment
    expect_identical(as.character(design(bcb)), c("~", "1"))
    expect_error(
        design(bcb) <- ~XXX,
        paste(
            "invalid class \"DESeqDataSet\" object:",
            "all variables in design formula must be columns in colData"
        )
    )
})
