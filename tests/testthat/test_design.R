context("design")

test_that("bcbioRNASeq", {
    # Not identical because of environment
    expect_equal(design(bcb), formula(~1))
    expect_error(
        design(bcb) <- formula(~XXX),
        paste(
            "invalid class \"DESeqDataSet\" object:",
            "all variables in design formula must be columns in colData"
        )
    )
})
