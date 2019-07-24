context("show")

test_that("bcbioRNASeq", {
    ## Stash fake metadata for code coverage.
    metadata(object)[["sampleMetadataFile"]] <- "XXX"
    metadata(object)[["gffFile"]] <- "XXX"
    output <- capture.output(object)
    ## Ensure that show method contains "bcbioRNASeq" in the first line.
    expect_true(
        grepl(
            pattern = "^bcbioRNASeq",
            x = output[[1L]]
        )
    )
})
